/// \file
/// \brief The crh_input_options_block class definitions

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
///
/// This program is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "crh_input_options_block.hpp"

#include <boost/algorithm/string/join.hpp>

#include "common/boost_addenda/program_options/layout_values_with_descs.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/program_options/prog_opt_num_range.hpp"

#include <limits>

using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::algorithm::join;
using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::numeric_limits;
using std::string;
using std::unique_ptr;

/// \brief The option name for the input file from which data should be read
const string crh_input_options_block::PO_INPUT_FILE_OR_STDIN    { "input-file-or-stdin"    };

/// \brief The option name for the format of the input data
const string crh_input_options_block::PO_INPUT_FORMAT           { "input-format"           };

/// \brief The option name for the minimum gap length to consider when parsing an alignment
const string crh_input_options_block::PO_MIN_GAP_LENGTH         { "min-gap-length"         };

/// \brief The option name for whether the code can assume that the input data is pre-grouped by query_id
const string crh_input_options_block::PO_INPUT_HITS_ARE_GROUPED { "input-hits-are-grouped" };

/// \brief A standard do_clone method
unique_ptr<options_block> crh_input_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string crh_input_options_block::do_get_block_name() const {
	return "Input";
}

/// \brief Add this block's options to the provided options_description
void crh_input_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                    ) {
	const auto &sep     = SUB_DESC_SEPARATOR;
	const auto &sub_sep = SUB_DESC_PAIR_SEPARATOR;

	const string format_varname { "<format>" };
	const string length_varname { "<length>" };

	const auto input_format_notifier           = [&] (const hits_input_format_tag &x) { the_spec.set_input_format          ( x ); };
	const auto min_gap_length_notifier         = [&] (const residx_t              &x) { the_spec.set_min_gap_length        ( x ); };
	const auto input_hits_are_grouped_notifier = [&] (const bool                  &x) { the_spec.set_input_hits_are_grouped( x ); };

	const str_vec input_format_descs = layout_values_with_descs(
		all_hits_input_format_tags,
		[] (const hits_input_format_tag &x) { return to_string( x ); },
		&description_of_input_format,
		sub_sep
	);

	arg_desc.add_options()
		(
			( PO_INPUT_FORMAT ).c_str(),
			value<hits_input_format_tag>()
				->value_name   ( format_varname                                 )
				->notifier     ( input_format_notifier                          )
				->default_value( crh_input_spec::DEFAULT_INPUT_FORMAT           ),
			( "Parse the input data from " + format_varname
				+ ", one of available formats:" + sep + join( input_format_descs, sep ) ).c_str()
		)
		(
			( PO_MIN_GAP_LENGTH ).c_str(),
			value< prog_opt_num_range<residx_t, 0, numeric_limits<residx_t>::max(), int64_t> >()
				->value_name   ( length_varname                                 )
				->notifier     ( min_gap_length_notifier                        )
				->default_value( crh_input_spec::DEFAULT_MIN_GAP_LENGTH         ),
			( "When parsing starts/stops from alignment data, ignore gaps of less than "
				+ length_varname + " residues").c_str()
		)
		(
			( PO_INPUT_HITS_ARE_GROUPED ).c_str(),
			bool_switch()
				->notifier     ( input_hits_are_grouped_notifier                )
				->default_value( crh_input_spec::DEFAULT_INPUT_HITS_ARE_GROUPED ),
			"Rely on the input hits being grouped by query protein"
			"\n(so the run is faster and uses less memory)"
		);

	static_assert( ! crh_input_spec::DEFAULT_READ_FROM_STDIN,        "If crh_input_spec::DEFAULT_READ_FROM_STDIN        isn't false, it might mess up the bool switch in here" );
	static_assert( ! crh_input_spec::DEFAULT_INPUT_HITS_ARE_GROUPED, "If crh_input_spec::DEFAULT_INPUT_HITS_ARE_GROUPED isn't false, it might mess up the bool switch in here" );
}

/// \brief Add a hidden option to the options_description for the input file
void crh_input_options_block::do_add_hidden_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                   ) {
	const string input_file_varname{ "<file>" };

	const auto input_file_or_stdin = [&] (const string &x) {
		if ( x == "-" ) {
			the_spec.set_read_from_stdin( true );
		}
		else {
			the_spec.set_input_file     ( x    );
		}
	};

	arg_desc.add_options()
		(
			PO_INPUT_FILE_OR_STDIN.c_str(),
			value<string>()
				->value_name   ( input_file_varname                                   )
				->notifier     ( input_file_or_stdin                                  ),
			( "Hidden option " + PO_INPUT_FILE_OR_STDIN + " with value " + input_file_varname + ".\n"
			  + "When " + input_file_varname + " is -, read from standard input." ).c_str()
		);
}

/// \brief Generate a description of any problem that makes the specified crh_input_options_block invalid
///        or none otherwise
str_opt crh_input_options_block::do_invalid_string(const variables_map &arg_variables_map ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                   ) const {

	if ( arg_variables_map.count( PO_MIN_GAP_LENGTH ) != 0 && ! arg_variables_map[ PO_MIN_GAP_LENGTH ].defaulted() ) {
		if ( the_spec.get_input_format() != hits_input_format_tag::HMMSEARCH_OUT ) {
			return "Cannot specify the minimum gap length for input formats that don't involve parsing gaps out of an alignment"s;
		}
	}
	return get_invalid_description( the_spec );
}

/// \brief Getter for the crh_input_spec that the crh_input_options_block configures
const crh_input_spec & crh_input_options_block::get_crh_input_spec() const {
	return the_spec;
}
