/// \file
/// \brief The crh_single_output_options_block class definitions

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

#include "crh_single_output_options_block.hpp"

#include <filesystem>

#include "cath/common/clone/make_uptr_clone.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::cath::rslv;

using ::boost::program_options::bool_switch;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::filesystem::path;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<options_block> crh_single_output_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string crh_single_output_options_block::do_get_block_name() const {
	return "Output";
}

/// \brief Add this block's options to the provided options_description
void crh_single_output_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                            const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                            ) {
	const string file_varname   { "<file>"   };

	const auto output_file_notifier          = [&] (const path &x) { the_spec.set_output_file          (           x ); };
	const auto summarise_notifier            = [&] (const bool &x) { the_spec.set_summarise            (           x ); };
	const auto generate_html_output_notifier = [&] (const bool &x) { the_spec.set_generate_html_output (           x ); };
	const auto json_output_notifier          = [&] (const bool &x) { the_spec.set_json_output          (           x ); };

	prm_desc.add_options()
		(
			string( PO_OUTPUT_FILE ).c_str(),
			value<path>()
				->value_name   ( file_varname                                         )
				->notifier     ( output_file_notifier                                 ),
			( "Write output to file " + file_varname + " (or, if unspecified, to stdout)" ).c_str()
		)
		(
			string( PO_SUMMARISE ).c_str(),
			bool_switch()
				->notifier     ( summarise_notifier                                   )
				->default_value( crh_single_output_spec::DEFAULT_SUMMARISE            ),
			"Output a brief text summary of the input data (rather than processing it)"
		)
		(
			string( PO_GENERATE_HTML_OUTPUT ).c_str(),
			bool_switch()
				->notifier     ( generate_html_output_notifier                        )
				->default_value( crh_single_output_spec::DEFAULT_GENERATE_HTML_OUTPUT ),
			"Output the results as HTML"
		)
		(
			string( PO_JSON_OUTPUT ).c_str(),
			bool_switch()
				->notifier     ( json_output_notifier                                 )
				->default_value( crh_single_output_spec::DEFAULT_JSON_OUTPUT          ),
			"Output the results as JSON"
		);

	static_assert( ! crh_single_output_spec::DEFAULT_GENERATE_HTML_OUTPUT,
		"If crh_single_output_spec::DEFAULT_GENERATE_HTML_OUTPUT isn't false, it might mess up the bool switch in here" );
	static_assert( ! crh_single_output_spec::DEFAULT_JSON_OUTPUT,
		"If crh_single_output_spec::DEFAULT_JSON_OUTPUT          isn't false, it might mess up the bool switch in here" );
}

/// \brief Generate a description of any problem that makes the specified crh_single_output_options_block invalid
///        or nullopt otherwise
str_opt crh_single_output_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                           ) const {
	return get_invalid_description( the_spec );
}

/// \brief Return all options names for this block
str_view_vec crh_single_output_options_block::do_get_all_options_names() const {
	return {
		PO_OUTPUT_FILE,
		PO_SUMMARISE,
		PO_GENERATE_HTML_OUTPUT,
		PO_JSON_OUTPUT,
	};
}

/// \brief Getter for the crh_single_output_spec that the crh_single_output_options_block configures
const crh_single_output_spec & crh_single_output_options_block::get_crh_single_output_spec() const {
	return the_spec;
}

/// \brief Get the crh_out_format for the crh_single_output_spec of the specified crh_single_output_options_block
///
/// \relates crh_single_output_options_block
crh_out_format cath::rslv::get_out_format(const crh_single_output_options_block &prm_crh_single_output_options_block ///< The crh_single_output_options_block to query
                                          ) {
	return get_out_format( prm_crh_single_output_options_block.get_crh_single_output_spec() );
}
