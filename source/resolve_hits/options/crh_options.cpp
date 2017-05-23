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

#include "common/boost_addenda/program_options/variables_map_contains.hpp"
#include "resolve_hits/options/spec/crh_spec.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;

using boost::none;
using boost::program_options::positional_options_description;
using boost::program_options::variables_map;
using std::map;
using std::pair;
using std::string;
using std::vector;

/// The name of the program that uses this executable_options
const string crh_options::PROGRAM_NAME("cath-resolve-hits");

/// \brief Get the name of the program that uses this executable_options
string crh_options::do_get_program_name() const {
	return PROGRAM_NAME;
}

/// \brief Get the positional options, which in this case is the input block's PO_INPUT_FILE_OR_STDIN option
positional_options_description crh_options::get_positional_options() {
	positional_options_description positionals;
	positionals.add( crh_input_options_block::PO_INPUT_FILE_OR_STDIN.c_str(), 1 );
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
	const auto &the_out_spec = the_single_output_ob.get_crh_single_output_spec();

	// If the user has specified neither an input file nor to read from stdin, then return an blank error string
	// (so the error will just be the basic "See 'cath-resolve-hits --help' for usage." message)
	if ( ! the_in_spec.get_input_file() && ! the_in_spec.get_read_from_stdin() && ! the_out_spec.get_export_css_file() ) {
		return string{};
	}

	const variables_map &vm           = get_variables_map();
	const auto          &input_format = get_crh_input_spec().get_input_format();
	if ( contains( vm, crh_score_options_block::PO_APPLY_CATH_RULES ) && ! vm[ crh_score_options_block::PO_APPLY_CATH_RULES ].defaulted() ) {
		if ( input_format != hits_input_format_tag::HMMER_DOMTBLOUT && input_format != hits_input_format_tag::HMMSEARCH_OUT ) {
			return "The --"
				+ crh_score_options_block::PO_APPLY_CATH_RULES
				+ " option cannot be used with the input format "
				+ to_string( input_format )
				+ "; CATH-Gene3D rules are only applied for "
				+ to_string( hits_input_format_tag::HMMER_DOMTBLOUT )
				+ " and "
				+ to_string( hits_input_format_tag::HMMSEARCH_OUT )
				+ " formats.";
		}
	}

	// Check that if the output_hmmsearch_aln option's enabled, the input format is HMMSEARCH_OUT
	if ( the_out_spec.get_output_hmmsearch_aln() && input_format != hits_input_format_tag::HMMSEARCH_OUT ) {
		return "Cannot use the --"
			+ crh_single_output_options_block::PO_OUTPUT_HMMSEARCH_ALN
			+ " option if using "
			+ to_string( input_format )
			+ " input format, must be using "
			+ to_string( hits_input_format_tag::HMMSEARCH_OUT );
	}

	// Store a map from score type to the equivalent "--worst-permissible-[...]" option name
	const auto worst_perm_opt_name_of_score = map<hit_score_type, string>{
		{ hit_score_type::FULL_EVALUE, crh_filter_options_block::PO_WORST_PERMISSIBLE_EVALUE   },
		{ hit_score_type::BITSCORE,    crh_filter_options_block::PO_WORST_PERMISSIBLE_BITSCORE },
		{ hit_score_type::CRH_SCORE,   crh_filter_options_block::PO_WORST_PERMISSIBLE_SCORE    },
	};

	// Store a map from the score type to the valid formats for which that "--worst-permissible-[...]" option may be specified
	const auto formats_for_worst_perm_opt_of_score = map< hit_score_type, hits_input_format_tag_vec >{
		{ hit_score_type::FULL_EVALUE, { hits_input_format_tag::RAW_WITH_EVALUES,                                      }, },
		{ hit_score_type::BITSCORE,    { hits_input_format_tag::HMMER_DOMTBLOUT, hits_input_format_tag::HMMSEARCH_OUT, }, },
		{ hit_score_type::CRH_SCORE,   { hits_input_format_tag::RAW_WITH_SCORES                                        }, },
	};
	//

	// For each such score type, check whether the relevant option has been specified and 
	// if so, validate that the input format is suitable or return an error string
	for (const auto &format_worse_conf : formats_for_worst_perm_opt_of_score) {
		const hit_score_type &score_type    = format_worse_conf.first;
		const auto           &valid_formats = format_worse_conf.second;
		const string         &option_name   = worst_perm_opt_name_of_score.at( score_type );
		
		
		if ( contains( vm, option_name ) && ! vm[ option_name ].defaulted() ) {
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

	if ( specifies_options_from_block( vm, crh_html_options_block{} ) && ( get_out_format( the_single_output_ob ) != crh_out_format::HTML ) ) {
		return
			"Cannot specify HTML options without setting the output format to HTML (with --"
			+ crh_single_output_options_block::PO_GENERATE_HTML_OUTPUT
			+ ")";
	}

	// If no error or help string, then return none
	return none;
}

/// \brief Get a string to prepend to the standard help
string crh_options::do_get_help_prefix_string() const {
	return "Usage: " + PROGRAM_NAME + R"( [options] <input_file>

)" + get_overview_string() + R"(

When <input_file> is -, the input is read from standard input.

The input data may contain unsorted hits for different query protein sequences.

However, if your input data is already grouped by query protein sequence, then
specify the --)" + crh_input_options_block::PO_INPUT_HITS_ARE_GROUPED + R"( flag for faster runs that use less memory.)";
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
				"Show help on the rules activated by the --"
					+ crh_score_options_block::PO_APPLY_CATH_RULES
					+ " option",
				get_crh_cath_rules_help_string()
			}
		},
	};
}

/// \brief Ctor, which initialises the detail_help_ob and adds the options_blocks to the parent executable_options
crh_options::crh_options() : detail_help_ob{ detail_help_spec() } {
	super::add_options_block( the_input_ob         );
	super::add_options_block( the_segment_ob       );
	super::add_options_block( the_score_ob         );
	super::add_options_block( the_filter_ob        );
	super::add_options_block( the_single_output_ob );
	super::add_options_block( the_html_ob          );
	super::add_options_block( detail_help_ob       );
}

/// \brief Build a crh_spec of the individual specs
crh_spec crh_options::get_crh_spec() const {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return crh_spec{
		the_input_ob.get_crh_input_spec(),
		the_segment_ob.get_crh_segment_spec(),
		the_score_ob.get_crh_score_spec(),
		the_filter_ob.get_crh_filter_spec(),
		the_single_output_ob.get_crh_single_output_spec(),
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
const crh_single_output_spec & crh_options::get_crh_single_output_spec() const {
	return the_single_output_ob.get_crh_single_output_spec();
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
	+ crh_score_options_block::PO_APPLY_CATH_RULES
	+ R"(]
------------------------------------

The --)"
	+ crh_score_options_block::PO_APPLY_CATH_RULES
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
)";
}