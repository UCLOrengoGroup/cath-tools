/// \file
/// \brief The cath_ssap_options class definitions

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
#include "cath_ssap_options.hpp"

#include <boost/program_options.hpp>
#include <boost/range/join.hpp>
#include <boost/shared_array.hpp>

#include "cath/acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "cath/acquirer/pdbs_acquirer/file_list_pdbs_acquirer.hpp"
#include "cath/acquirer/pdbs_acquirer/istream_pdbs_acquirer.hpp"
#include "cath/acquirer/selection_policy_acquirer/selection_policy_acquirer.hpp"
#include "cath/alignment/alignment.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.hpp"
#include "cath/chopping/domain/domain.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/argc_argv_faker.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/options/options_block/options_block.hpp"
#include "cath/outputter/superposition_outputter/superposition_outputter.hpp"
#include "cath/outputter/superposition_outputter/superposition_outputter_list.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_source_file_set.hpp"
#include "cath/superposition/superposition_context.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;

using boost::filesystem::path;
using boost::none;
using boost::program_options::positional_options_description;
using boost::range::join;
using std::string;

/// \brief TODOCUMENT
const string cath_ssap_options::PO_CITATION_HELP{ "citation-help" };

/// \brief The name of the program that uses this executable_options
const string cath_ssap_options::PROGRAM_NAME    { "cath-ssap"     };

/// \brief Get the options for the "Detailed Help" block
str_str_str_pair_map cath_ssap_options::detail_help_spec() {
	return {
		{ "alignment-help",         { "Help on alignment format",                 get_ssap_alignment_format_help_string() } },
		{ "scores-help",            { "Help on scores format",                    get_ssap_matches_format_help_string()   } },
		{ PO_CITATION_HELP.c_str(), { "Help on SSAP authorship & how to cite it", get_ssap_citation_help_string()         } }
	};
}

/// \brief Get the name of the program that uses this executable_options
string cath_ssap_options::do_get_program_name() const {
	return PROGRAM_NAME;
}

/// \brief Specify the name option in old_ssap_options_block as the positional option
///
/// This overrides a virtual method in executable_options
positional_options_description cath_ssap_options::get_positional_options() {
	positional_options_description positionals;
	positionals.add( old_ssap_options_block::PO_NAME.c_str(), -1 );
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
str_opt cath_ssap_options::do_get_error_or_help_string() const {
	// If detailed help was requested, then provide it
	if ( the_detail_help_options_block.has_help_string() ) {
		return the_detail_help_options_block.help_string();
	}

	// If there are no proteins were specified, just output the standard usage error string
	if ( ! get_old_ssap_options().protein_names_specified() ) {
		return ""s;
	}

	// If any of the required files are not valid input files then complain
	const bool    has_superposition   = has_superposition_dir( the_ssap_options_block );

	const auto    protein_file_types  = the_ssap_options_block.get_protein_source_files()->get_file_set();
	const auto    supn_file_type      = has_superposition ? data_file_vec{} : data_file_vec{ 1, data_file::PDB };
	const auto    required_file_types = copy_build<data_file_set>( join( protein_file_types, supn_file_type ) );
	const str_vec protein_names       = { the_ssap_options_block.get_protein_name_a(),
	                                      the_ssap_options_block.get_protein_name_b() };
	const auto    &the_data_dirs_spec = get_data_dirs_spec();
	path_vec required_input_files;
	try {
		for (const data_file &required_file_type : required_file_types) {
			for (const string &protein_name : protein_names) {
				required_input_files.push_back( find_file( the_data_dirs_spec, required_file_type, protein_name ) );
			}
		}
	}
	catch (const runtime_error_exception &the_exception) {
		return "Problem with input file: "s + the_exception.what();
	}
	for (const path &required_input_file : required_input_files) {
		if ( ! options_block::is_acceptable_input_file( required_input_file ) ) {
			return "Problem with input file: Required input file "
			       + required_input_file.string()
			       + " is not a valid, non-empty input file";
		}
	}

	if ( get_domains().size() > 2 ) {
		return "Cannot specify regions more than twice for cath-ssap"s;
	}

	return none;
}

/// \brief Get a string to prepend to the standard help
string cath_ssap_options::do_get_help_prefix_string() const {
	return "Usage: " + PROGRAM_NAME + R"( [options] <protein1> <protein2>

)" + get_overview_string() + R"(

)" + PROGRAM_NAME + R"( uses two types of structural comparison:
  1. Fast SSAP: a quick secondary-structure based SSAP alignment
  2. Slow SSAP: residue alignment only

If both structures have more than one SS element, a fast SSAP is run first.)"
	" If the fast SSAP score isn't good, another fast SSAP is run with looser cutoffs."
	" If the (best) fast SSAP score isn't good, a slow SSAP is run."
	" Only the best of these scores is output."
	" These behaviours can be configured using the parameters below.)";
}

/// \brief Get a string to append to the standard help (just empty here)
string cath_ssap_options::do_get_help_suffix_string() const {
	return "";
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string cath_ssap_options::do_get_overview_string() const {
	return R"(Run a SSAP pairwise structural alignment
[algorithm devised by C A Orengo and W R Taylor, see --)" + PO_CITATION_HELP + "]";
}

/// \brief Check that these options are OK to use
///
/// \pre This should only be called when they are OK to use,
///      else this will throw an invalid_argument_exception
void cath_ssap_options::check_ok_to_use() const {
	if ( get_error_or_help_string() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use invalid cath_ssap_options"));
	}
}

/// \brief Default ctor for cath_ssap_options
///
/// This adds the options blocks to the parent executable_options class
cath_ssap_options::cath_ssap_options() : the_detail_help_options_block( detail_help_spec() ) {
	super::add_options_block( the_ssap_options_block        );
	super::add_options_block( the_data_dirs_options_block   );
	super::add_options_block( the_align_regions_ob          );
	super::add_options_block( the_detail_help_options_block );
}

/// \brief A getter for the old_ssap_options_block
const old_ssap_options_block & cath_ssap_options::get_old_ssap_options() const {
	return the_ssap_options_block;
}

/// \brief A getter for the data_dirs_options_block
const data_dirs_spec & cath_ssap_options::get_data_dirs_spec() const {
	return the_data_dirs_options_block.get_data_dirs_spec();
}

/// \brief A getter for the cath_ssap_options
const domain_vec & cath_ssap_options::get_domains() const {
	return the_align_regions_ob.get_align_domains();
}

/// \brief Return the string containing help on the SSAP matches format
string cath::opts::get_ssap_matches_format_help_string() {
	return R"(Help on Standard SSAP Scores Output
===================================

Overview
--------
  Format  : "%6s  %6s %4d %4d %6.2f %4d %4d %4d %6.2f"
  Columns : prot1 prot2 length1 length2 ssap_score num_equivs overlap_pc seq_id_pc rmsd

Column descriptions
------------------
  Column 1 (prot1     ) : Name of protein 1
  Column 2 (prot2     ) : Name of protein 2
  Column 3 (length1   ) : Length of protein 1
  Column 4 (length2   ) : Length of protein 2
  Column 5 (ssap_score) : SSAP score for structural comparison (0-100)
  Column 6 (num_equivs) : Number of equivalent/aligned residues
  Column 7 (overlap_pc) : Percentage overlap  (100% x overlap /length of largest)
  Column 8 (seq_id_pc ) : Percentage identity (100% x identity/length of smallest)
  Column 9 (rmsd      ) : RMSD of superposed structures

Example
-------
  "1cukA03  1hjpA03   48   44  94.92   44   91   97   0.71"
)";
}

/// \brief Return the string containing help on the SSAP alignment format
string cath::opts::get_ssap_alignment_format_help_string() {
	return R"(Output Format: Standard SSAP Alignment Format
=============================================

Format
------
String Format "%4d %c %c %c  %3d  %c %c %c %4d"

Column 1: Protein 1 PDB residue number (excluding insert character)
Column 2: Protein 1 Secondary structure character
Column 3: Protein 1 PDB residue insert character
Column 4: Protein 1 Residue code (one letter code)
Column 5: SSAP residue score (0-100)
Column 6: Protein 2 Residue code (one letter code)
Column 7: Protein 2 PDB residue insert character
Column 8: Protein 2 Secondary structure character
Column 9: Protein 2 PDB residue number (excluding insert character)

Discussion
----------
Alignment files are named "prot1prot2.list" e.g. 1cukA031hjp03.list

Secondary structure characters are taken from the wolf file (modified DSSP)
These are a subset of the DSSP definitions (Kabsch and Sander, 1983)

 * E - extended strand, participates in beta-ladder
 * G - 3-helix (3/10 helix)
 * H - 4-helix (alpha-helix)
 * T - hydrogen-bonded turn

Example
-------
 158 H 0 D    0  0 0 0    0
 159 H 0 A   85  A 0 0  159
 160 H 0 E   88  E 0 H  160
 161 H 0 Q   89  Q 0 H  161
 162 H 0 E   90  E 0 H  162
 163 H 0 A   93  A 0 H  163
 164 H 0 V   92  V 0 H  164
)";
}

/// \brief Return the string containing help on the SSAP alignment format
string cath::opts::get_ssap_citation_help_string() {
	return R"(Please cite: "Protein Structure Alignment", Taylor and Orengo [1989]
Journal of Molecular Biology 208, 1-22
PMID: 2769748

Many people have contributed to this code, most notably:
  * Tony E Lewis               (  2011 - ....)
  * Oliver C Redfern           (~ 2003 - 2011)
  * James E Bray, Ian Sillitoe (~ 2000 - 2003)
  * Andrew C R Martin          (considerable edits around 2001)
)";
}