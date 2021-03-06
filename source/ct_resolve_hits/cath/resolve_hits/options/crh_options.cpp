/// \file
/// \brief The crh_options class definitions

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

#include "crh_options.hpp"

#include <fmt/core.h>

#include "cath/common/boost_addenda/program_options/variables_map_contains.hpp"
#include "cath/resolve_hits/options/spec/crh_spec.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::cath::rslv;

using ::boost::program_options::positional_options_description;
using ::boost::program_options::variables_map;
using ::std::map;
using ::std::nullopt;
using ::std::string;
using ::std::string_view;

/// \brief Get the name of the program that uses this executable_options
string_view crh_options::do_get_program_name() const {
	return PROGRAM_NAME;
}

/// \brief Get the positional options, which in this case is the input block's PO_INPUT_FILE_OR_STDIN option
positional_options_description crh_options::get_positional_options() {
	positional_options_description positionals;
	positionals.add( string( crh_input_options_block::PO_INPUT_FILE_OR_STDIN ).c_str(), 1 );
	return positionals;
}

/// \brief Review all specified options and return a string containing any errors or a help string
///        (possibly using a description of all visible options)
///
/// This is a concrete definition of a virtual method that's pure in executable_options
///
/// This should only be called by executable_options, as the last step of the parse_options()
/// method, after all real parsing has completed.
///
/// \pre The options should have been parsed
///
/// \returns Any error/help string arising from the newly specified options
///          or an empty string if there aren't any
str_opt crh_options::do_get_error_or_help_string() const {
	// If detailed help was requested, then provide it
	if ( detail_help_ob.has_help_string() ) {
		return detail_help_ob.help_string();
	}

	const auto &the_in_spec  = the_input_ob.get_crh_input_spec();
	const auto &the_out_spec = the_output_ob.get_crh_output_spec();

	// If the user has specified neither an input file nor to read from stdin, then return an blank error string
	// (so the error will just be the basic "See 'cath-resolve-hits --help' for usage." message)
	if ( ! the_in_spec.get_input_file() && ! the_in_spec.get_read_from_stdin() && ! the_out_spec.get_export_css_file() ) {
		return string{};
	}

	const variables_map &local_vm           = get_variables_map();
	const auto          &input_format = get_crh_input_spec().get_input_format();
	if ( specifies_option( local_vm, string( crh_score_options_block::PO_APPLY_CATH_RULES ) ) ) {
		if ( input_format != hits_input_format_tag::HMMER_DOMTBLOUT && input_format != hits_input_format_tag::HMMSEARCH_OUT ) {
			return ::fmt::format( "The --{} option cannot be used with the input format {}; CATH-Gene3D rules are only "
			                      "applied for {} and {} formats.",
			                      crh_score_options_block::PO_APPLY_CATH_RULES,
			                      to_string( input_format ),
			                      to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
			                      to_string( hits_input_format_tag::HMMSEARCH_OUT ) );
		}
	}

	// Check that if the output_hmmer_aln option's enabled, the input format is HMMSCAN_OUT / HMMSEARCH_OUT
	if ( the_out_spec.get_output_hmmer_aln() && input_format != hits_input_format_tag::HMMSCAN_OUT  && input_format != hits_input_format_tag::HMMSEARCH_OUT ) {
		return ::fmt::format( "Cannot use the --{} option if using {} input format, must be using {} or ",
		                      crh_output_options_block::PO_OUTPUT_HMMER_ALN,
		                      to_string( input_format ),
		                      to_string( hits_input_format_tag::HMMSCAN_OUT ),
		                      to_string( hits_input_format_tag::HMMSEARCH_OUT ) );
	}

	// Store a map from score type to the equivalent "--worst-permissible-[...]" option name
	const auto worst_perm_opt_name_of_score = map<hit_score_type, string_view>{
		{ hit_score_type::FULL_EVALUE, crh_filter_options_block::PO_WORST_PERMISSIBLE_EVALUE   },
		{ hit_score_type::BITSCORE,    crh_filter_options_block::PO_WORST_PERMISSIBLE_BITSCORE },
		{ hit_score_type::CRH_SCORE,   crh_filter_options_block::PO_WORST_PERMISSIBLE_SCORE    },
	};

	// Store a map from the score type to the valid formats for which that "--worst-permissible-[...]" option may be specified
	const auto formats_for_worst_perm_opt_of_score = map< hit_score_type, hits_input_format_tag_vec >{
		{ hit_score_type::FULL_EVALUE, { hits_input_format_tag::RAW_WITH_EVALUES,                                                                          }, },
		{ hit_score_type::BITSCORE,    { hits_input_format_tag::HMMER_DOMTBLOUT, hits_input_format_tag::HMMSCAN_OUT, hits_input_format_tag::HMMSEARCH_OUT, }, },
		{ hit_score_type::CRH_SCORE,   { hits_input_format_tag::RAW_WITH_SCORES                                                                            }, },
	};
	//

	// For each such score type, check whether the relevant option has been specified and
	// if so, validate that the input format is suitable or return an error string
	for (const auto &format_worse_conf : formats_for_worst_perm_opt_of_score) {
		const hit_score_type &score_type    = format_worse_conf.first;
		const auto           &valid_formats = format_worse_conf.second;
		const string_view    &option_name   = worst_perm_opt_name_of_score.at( score_type );

		if ( specifies_option( local_vm, string( option_name ) ) ) {
			if ( ! contains( valid_formats, input_format ) ) {
				return "Cannot set worst permissible "
					+ to_string( score_type   )
					+ " for input format "
					+ to_string( input_format )
					+ ", for which "
					+ to_string( score_type   )
					+ " isn't the primary score type.";
			}
		}
	}

	if ( specifies_options_from_block( local_vm, crh_html_options_block{} ) && ! has_any_html_output( the_output_ob ) ) {
		return ::fmt::format( "Cannot specify HTML options without outputting any HTML (with --{})",
		                      crh_output_options_block::PO_HTML_OUTPUT_TO_FILE );
	}

	// Complain if using MIN_HMM_COVERAGE / MIN_DC_HMM_COVERAGE for any format other than HMMSEARCH_OUT
	if ( input_format != hits_input_format_tag::HMMSEARCH_OUT ) {
		for (const string_view &hmm_cov_opt : { crh_filter_options_block::PO_MIN_HMM_COVERAGE,
		                                        crh_filter_options_block::PO_MIN_DC_HMM_COVERAGE } ) {
			const string hmm_cov_opt_str( hmm_cov_opt );
			if ( specifies_option( local_vm, hmm_cov_opt_str ) ) {
				return "Cannot specify --"
					+ hmm_cov_opt_str
					+ " unless using input format "
					+ to_string( hits_input_format_tag::HMMSEARCH_OUT );
			}
		}
	}

	// If no error or help string, then return nullopt
	return nullopt;
}

/// \brief Get a string to prepend to the standard help
string crh_options::do_get_help_prefix_string() const {
	return ::fmt::format(
	  R"(Usage: {} [options] <input_file>

{}

When <input_file> is -, the input is read from standard input.

The input data may contain unsorted hits for different query protein sequences.

However, if your input data is already grouped by query protein sequence, then
specify the --{} flag for faster runs that use less memory.)",
	  PROGRAM_NAME,
	  get_overview_string(),
	  crh_input_options_block::PO_INPUT_HITS_ARE_GROUPED );
}

/// \brief Get a string to append to the standard help
string crh_options::do_get_help_suffix_string() const {
	return "\nThe standard output is one line per selected hit,"
		" preceded by header lines (beginning \"#\"), the last of which (beginning \"#FIELDS\") lists the fields in the file, typically:\n"
		"  #FIELDS query-id match-id score boundaries resolved\n"
		"(`boundaries` and `resolved` describe a domain's starts / stops; `resolved` may include adjustments made to resolve overlaps between hits)\n";
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string crh_options::do_get_overview_string() const {
	return R"(Collapse a list of domain matches to your query sequence(s) down to the
non-overlapping subset (ie domain architecture) that maximises the sum of the
hits' scores.)";
}

/// \brief Get the options for the "Detailed Help" block
str_str_str_pair_map crh_options::detail_help_spec() {
	return {
		{
			"raw-format-help",
			{
				"Show help about the raw input formats ("
					+ to_string( hits_input_format_tag::RAW_WITH_SCORES )
					+ " and "
					+ to_string( hits_input_format_tag::RAW_WITH_EVALUES )
					+ ")",
				get_crh_raw_format_help_string()
			}
		},
		{
			"cath-rules-help",
			{
				::fmt::format(
					"Show help on the rules activated by the (DEPRECATED) --{} option",
					crh_score_options_block::PO_APPLY_CATH_RULES
				),
				get_crh_cath_rules_help_string()
			}
		},
	};
}

/// \brief Ctor, which initialises the detail_help_ob and adds the options_blocks to the parent executable_options
crh_options::crh_options() : detail_help_ob{ detail_help_spec() } {
	super::add_options_block( the_input_ob   );
	super::add_options_block( the_segment_ob );
	super::add_options_block( the_score_ob   );
	super::add_options_block( the_filter_ob  );
	super::add_options_block( the_output_ob  );
	super::add_options_block( the_html_ob    );
	super::add_options_block( detail_help_ob );
}

/// \brief Build a crh_spec of the individual specs
crh_spec crh_options::get_crh_spec() const {
	return crh_spec{
		the_input_ob.get_crh_input_spec(),
		the_segment_ob.get_crh_segment_spec(),
		the_score_ob.get_crh_score_spec(),
		the_filter_ob.get_crh_filter_spec(),
		get_deprecated_single_output_spec( the_output_ob ),
		the_output_ob.get_crh_output_spec(),
		the_html_ob.get_crh_html_spec()
	};
}

/// \brief Getter for the cath-resolve-hits input options_block
const crh_input_spec & crh_options::get_crh_input_spec() const {
	return the_input_ob.get_crh_input_spec();
}

/// \brief Getter for the cath-resolve-hits segment options_block
const crh_segment_spec & crh_options::get_crh_segment_spec() const {
	return the_segment_ob.get_crh_segment_spec();
}

/// \brief Getter for the cath-resolve-hits score options_block
const crh_score_spec & crh_options::get_crh_score_spec() const {
	return the_score_ob.get_crh_score_spec();
}

/// \brief Getter for the cath-resolve-hits filter options_block
const crh_filter_spec & crh_options::get_crh_filter_spec() const {
	return the_filter_ob.get_crh_filter_spec();
}

/// \brief Getter for the cath-resolve-hits output options_block
const crh_output_spec & crh_options::get_crh_output_spec() const {
	return the_output_ob.get_crh_output_spec();
}

/// \brief Getter for the cath-resolve-hits html options_block
const crh_html_spec & crh_options::get_crh_html_spec() const {
	return the_html_ob.get_crh_html_spec();
}

/// \brief Get the text for the raw format help
string cath::rslv::get_crh_raw_format_help_string() {
	return  R"(Format for "raw" input data
---------------------------

One hit per line, using the following space-separated fields:

 1. query_protein_id : An identifier for the query protein sequence
 2. match_id         : An identifier for the match
 3. score            : A (strictly positive) score indicating how good it would be to have that hit in the final result
 4. starts_stops     : The starts/stops on the query sequence, given in the format: 37-124,239-331

Example lines:

qyikaz 1mkfA01/12-210-i5_1,2.9e-20 2983.29780221221 3-103
qyikaz 1mkfA01/12-210-i5_2,4.9e-19 3510.41568607646 101-224
qyikaz 1mkfA01/12-210-i5_3,7e-25 3552.10980383852 825-928
qyikaz 1mkfA01/12-210-i5_4,3.5e-15 2470.04912752062 953-1053
)";
}

/// \brief Get the text for the apply CATH-Gene3d rules help
string cath::rslv::get_crh_cath_rules_help_string() {
	return R"(CATH Rules Help [--)"
	+ string( crh_score_options_block::PO_APPLY_CATH_RULES )
	+ R"(]
------------------------------------

[*DEPRECATED*] (please raise a GitHub issue if you want --)"
	+ string( crh_score_options_block::PO_APPLY_CATH_RULES )
	+ R"( to be kept)

The --)"
	+ string( crh_score_options_block::PO_APPLY_CATH_RULES )
	+ R"( option applies the following CATH-Gene3D specific rules when parsing from )"
	+ to_string( hits_input_format_tag::HMMER_DOMTBLOUT )
	+ R"( or )"
	+ to_string( hits_input_format_tag::HMMSEARCH_OUT )
	+ R"( files.

If hit's match ID is like "dc_72a964d791dea7a3dd35a8bbf49385b8" (matches /^dc_\w{32}$/):
 * use the ali_from/ali_to fields rather than env_from/env_to to determine the final start/stop
 * ignore gaps when parsing an alignment from a )"
	+ to_string( hits_input_format_tag::HMMSEARCH_OUT )
	+ R"(file (ie keep the hit as one continuous segment)

If the conditional-evalue is <= 0.001 but the independent-value is > 0.001, then quarter the bitscore when parsing the hit.

[*DEPRECATED*] (please raise a GitHub issue if you want --)"
	+ string( crh_score_options_block::PO_APPLY_CATH_RULES )
	+ R"( to be kept)
)";
}