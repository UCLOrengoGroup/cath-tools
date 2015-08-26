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
#include "cath_ssap_options.h"

#include <boost/program_options.hpp>
#include <boost/range/join.hpp>
#include <boost/shared_array.hpp>

#include "alignment/alignment.h"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.h"
#include "common/algorithm/copy_build.h"
#include "common/argc_argv_faker.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"
#include "exception/not_implemented_exception.h"
#include "exception/runtime_error_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "options/acquirer/alignment_acquirer/alignment_acquirer.h"
#include "options/acquirer/pdbs_acquirer/file_list_pdbs_acquirer.h"
#include "options/acquirer/pdbs_acquirer/istream_pdbs_acquirer.h"
#include "options/acquirer/selection_policy_acquirer/selection_policy_acquirer.h"
#include "options/acquirer/superposition_acquirer/align_based_superposition_acquirer.h"
#include "options/options_block/options_block.h"
#include "options/outputter/superposition_outputter/superposition_outputter.h"
#include "options/outputter/superposition_outputter/superposition_outputter_list.h"
#include "structure/protein/protein_source_file_set/protein_source_file_set.h"
#include "superposition/superposition_context.h"

#include <iostream>

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace std;

using boost::range::join;

/// \brief TODOCUMENT
const str_str_str_pair_map cath_ssap_options::DETAIL_HELP_SPEC = {
	{ "alignment-help", { "Help on alignment format", get_ssap_alignment_format_help_string() } },
	{ "scores-help",    { "Help on scores format",    get_ssap_matches_format_help_string()   } }
};

/// \brief TODOCUMENT
const string cath_ssap_options::PROGRAM_NAME( "cath-ssap" );

/// \brief Provide a name for the executable
///
/// This is a concrete definition of a virtual method that's pure in executable_options
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
string cath_ssap_options::do_update_error_or_help_string(const options_description &arg_visible_program_options ///< The full options_description of visible options
                                                         ) const {
	// If help was requested, then provide it
	if (the_misc_options_block.get_help()) {
		return the_misc_options_block.get_help_string(arg_visible_program_options, get_usage_string(), "");
	}

	// If version information was requested, then provide it
	if (the_misc_options_block.get_version()) {
		return the_misc_options_block.get_version_string(get_program_name(), get_version_description_string());
	}

	// If detailed help was requested, then provide it
	if (the_detail_help_options_block.has_help_string()) {
		return the_detail_help_options_block.help_string();
	}

	// If there are no proteins were specified, just output the standard usage error string
	if ( ! get_old_ssap_options().protein_names_specified() ) {
		return get_standard_usage_error_string();
	}

	// If any of the required files are not valid input files then complain
	const bool    has_superposition   = has_superposition_dir( the_ssap_options_block );

	const auto    protein_file_types  = the_ssap_options_block.get_protein_source_files()->get_file_set();
	const auto    supn_file_type      = has_superposition ? data_file_vec{} : data_file_vec{ 1, data_file::PDB };
	const auto    required_file_types = copy_build<data_file_set>( join( protein_file_types, supn_file_type ) );
	const str_vec protein_names       = { the_ssap_options_block.get_protein_name_a(),
	                                      the_ssap_options_block.get_protein_name_b() };
	path_vec required_input_files;
	try {
		for (const data_file &required_file_type : required_file_types) {
			for (const string &protein_name : protein_names) {
				required_input_files.push_back( find_file( the_data_dirs_options_block, required_file_type, protein_name ) );
			}
		}
	}
	catch (const runtime_error_exception &the_exception) {
		return string( "Problem with input file: " ) + the_exception.what() + "\n" + get_standard_usage_error_string();
	}
	for (const path &required_input_file : required_input_files) {
		if ( ! options_block::is_acceptable_input_file( required_input_file ) ) {
			return "Problem with input file: Required input file "
			       + required_input_file.string()
			       + " is not a valid, non-empty input file\n"
			       + get_standard_usage_error_string();
		}
	}

	return "";
}

/// \brief Output the standard program overview to an ostream
string cath_ssap_options::get_version_description_string() {
	ostringstream version_ss;
	version_ss << R"(Protein structure comparison algorithm
Devised by Christine A Orengo and William R Taylor

Please cite: "Protein Structure Alignment", Taylor and Orengo [1989]
             Journal of Molecular Biology 208, 1-22
             PMID: 2769748

Many people have contributed to this code, most notably:
  * Tony E Lewis               (  2011 - ....)
  * Oliver C Redfern           (~ 2003 - 2011)
  * James E Bray, Ian Sillitoe (~ 2000 - 2003)
  * Andrew C R Martin          (considerable edits around 2001))";
	return version_ss.str();
}

/// \brief Return the standard usage string
string cath_ssap_options::get_usage_string() {
	ostringstream usage_ss;
	usage_ss
		<< "Usage: " << PROGRAM_NAME << R"( [options] <protein1> <protein2>

Run a SSAP pairwise structural alignment

)" << PROGRAM_NAME << R"( uses two types of structural comparison:
  1. Fast SSAP: a quick secondary-structure based SSAP alignment
  2. Slow SSAP: residue alignment only

If both structures have more than one SS element, a fast SSAP is run first. If the fast SSAP score isn't good, another fast SSAP is run with looser cutoffs. If the (best) fast SSAP score isn't good, a slow SSAP is run. Only the best of these scores is output. These behaviours can be configured using the parameters below.)";
	return usage_ss.str();
}

/// \brief Check that these options are OK to use
///
/// \pre This should only be called when they are OK to use,
///      else this will throw an invalid_argument_exception
void cath_ssap_options::check_ok_to_use() const {
	if (!get_error_or_help_string().empty()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use invalid cath_ssap_options"));
	}
	if (the_misc_options_block.get_help()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use cath_ssap_options that is set to do_nothing"));
	}
}

/// \brief Default ctor for cath_ssap_options
///
/// This adds the options blocks to the parent executable_options class
cath_ssap_options::cath_ssap_options() : the_detail_help_options_block(DETAIL_HELP_SPEC) {
	super::add_options_block( the_misc_options_block        );
	super::add_options_block( the_ssap_options_block        );
	super::add_options_block( the_data_dirs_options_block   );
	super::add_options_block( the_detail_help_options_block );
}

/// \brief A getter for the old_ssap_options_block
const old_ssap_options_block cath_ssap_options::get_old_ssap_options() const {
	return the_ssap_options_block;
}

/// \brief A getter for the data_dirs_options_block
const data_dirs_options_block cath_ssap_options::get_data_dirs_options() const {
	return the_data_dirs_options_block;
}

/// \brief Return the string containing help on the SSAP matches format
string cath::opts::get_ssap_matches_format_help_string() {
	ostringstream ssap_matches_format_help_ss;
	ssap_matches_format_help_ss << R"(Help on Standard SSAP Scores Output
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
  "1cukA03  1hjpA03   48   44  94.92   44   91   97   0.71")";
	return ssap_matches_format_help_ss.str();
}


/// \brief Return the string containing help on the SSAP alignment format
string cath::opts::get_ssap_alignment_format_help_string() {
	ostringstream ssap_matches_format_help_ss;
	ssap_matches_format_help_ss << R"(Output Format: Standard SSAP Alignment Format
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
 164 H 0 V   92  V 0 H  164)";

	return ssap_matches_format_help_ss.str();
}
