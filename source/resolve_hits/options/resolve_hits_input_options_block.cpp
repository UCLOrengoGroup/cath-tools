/// \file
/// \brief The resolve_hits_input_options_block class definitions

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

#include "resolve_hits_input_options_block.h"

#include "common/clone/make_uptr_clone.h"

using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using boost::program_options::options_description;
using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief TODOCUMENT
const string resolve_hits_input_options_block::PO_INPUT_FILE_OR_STDIN   { "input-file-or-stdin"       };

/// \brief TODOCUMENT
const string resolve_hits_input_options_block::PO_INPUT_HITS_ARE_GROUPED{ "assume-input-hits-grouped" };

/// \brief TODOCUMENT
const string resolve_hits_input_options_block::PO_DOMTBLOUT             { "domtblout-input"          };

/// \brief TODOCUMENT
const string resolve_hits_input_options_block::PO_HMMERALN              { "hmmeraln-input"           };

/// \brief TODOCUMENT
const string resolve_hits_input_options_block::PO_RAW_SCORE_IS_EVALUE   { "raw-score-is-evalue"      };

/// \brief TODOCUMENT
const string resolve_hits_input_options_block::PO_APPLY_CATH_RULES      { "apply-cath-rules"         };

/// \brief A standard do_clone method.
unique_ptr<options_block> resolve_hits_input_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string resolve_hits_input_options_block::do_get_block_name() const {
	return "Resolve hits";
}

/// \brief Add this block's options to the provided options_description
void resolve_hits_input_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                             ) {
	const auto input_hits_are_sorted_notifier = [&] (const bool &x) { the_spec.set_input_hits_are_grouped( x ); };
	const auto format_is_domtblout_notifier   = [&] (const bool &x) { the_spec.set_domtblout             ( x ); };
	const auto format_is_hmmeraln_notifier    = [&] (const bool &x) { the_spec.set_hmmeraln              ( x ); };
	const auto raw_score_is_evalue_notifier   = [&] (const bool &x) { the_spec.set_raw_score_is_evalue   ( x ); };
	const auto apply_cath_rules_notifier      = [&] (const bool &x) { the_spec.set_apply_cath_rules      ( x ); };

	arg_desc.add_options()
		(
			( PO_INPUT_HITS_ARE_GROUPED ).c_str(),
			bool_switch()
				->notifier     ( input_hits_are_sorted_notifier                          )
				->default_value( resolve_hits_input_spec::DEFAULT_INPUT_HITS_ARE_GROUPED ),
			"Assume that the input hits are grouped by query protein"
			" (so the run is faster and uses less memory)."
		)
		(
			( PO_DOMTBLOUT ).c_str(),
			bool_switch()
				->notifier     ( format_is_domtblout_notifier                         )
				->default_value( resolve_hits_input_spec::DEFAULT_DOMTBLOUT           ),
			"Assume the input is the HMMER domtblout fomat."
		)
		(
			( PO_HMMERALN ).c_str(),
			bool_switch()
				->notifier     ( format_is_hmmeraln_notifier                          )
				->default_value( resolve_hits_input_spec::DEFAULT_HMMERALN            ),
			"Assume the input is the hmmsearch output fomat."
			"."
		)
		(
			( PO_RAW_SCORE_IS_EVALUE ).c_str(),
			bool_switch()
				->notifier     ( raw_score_is_evalue_notifier                         )
				->default_value( resolve_hits_input_spec::DEFAULT_RAW_SCORE_IS_EVALUE ),
			"Read the scores from default \"raw\" format as evalues and convert to scores accordingly."
		)
		(
			( PO_APPLY_CATH_RULES ).c_str(),
			bool_switch()
				->notifier     ( apply_cath_rules_notifier                            )
				->default_value( resolve_hits_input_spec::DEFAULT_APPLY_CATH_RULES    ),
			"Apply CATH-Gene3D specific rules to the parsing and processing."
		);

	static_assert( ! resolve_hits_input_spec::DEFAULT_DOMTBLOUT,           "If resolve_hits_input_spec::DEFAULT_DOMTBLOUT           isn't false, it might mess up the bool switch in here" );
	static_assert( ! resolve_hits_input_spec::DEFAULT_HMMERALN,            "If resolve_hits_input_spec::DEFAULT_HMMERALN            isn't false, it might mess up the bool switch in here" );
	static_assert( ! resolve_hits_input_spec::DEFAULT_RAW_SCORE_IS_EVALUE, "If resolve_hits_input_spec::DEFAULT_RAW_SCORE_IS_EVALUE isn't false, it might mess up the bool switch in here" );
	static_assert( ! resolve_hits_input_spec::DEFAULT_APPLY_CATH_RULES,    "If resolve_hits_input_spec::DEFAULT_APPLY_CATH_RULES    isn't false, it might mess up the bool switch in here" );
}

/// \brief TODOCUMENT
void resolve_hits_input_options_block::do_add_hidden_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
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
			( "TODOCUMENT option " + PO_INPUT_FILE_OR_STDIN + " with value " + input_file_varname + ".\n"
			  + "When " + input_file_varname + " is -, read from standard input." ).c_str()
		);
}

/// \brief TODOCUMENT
opt_str resolve_hits_input_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                            ) const {

	return get_invalid_description( the_spec );
}

/// \brief TODOCUMENT
const resolve_hits_input_spec & resolve_hits_input_options_block::get_resolve_hits_input_spec() const {
	return the_spec;
}
