/// \file
/// \brief The cath_superpose_options class definitions

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

#include "cath_superpose_options.h"

#include <boost/program_options.hpp>

#include "alignment/alignment.h"
#include "alignment/common_atom_selection_policy/common_atom_select_ca_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.h"
#include "common/logger.h"
#include "common/type_aliases.h"
#include "display/display_spec/display_spec.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "options/acquirer/alignment_acquirer/alignment_acquirer.h"
#include "options/acquirer/pdbs_acquirer/pdbs_acquirer.h"
#include "options/acquirer/selection_policy_acquirer/selection_policy_acquirer.h"
#include "options/acquirer/superposition_acquirer/align_based_superposition_acquirer.h"
#include "options/outputter/alignment_outputter/alignment_outputter.h"
#include "options/outputter/alignment_outputter/alignment_outputter_list.h"
#include "options/outputter/superposition_outputter/superposition_outputter.h"
#include "options/outputter/superposition_outputter/superposition_outputter_list.h"
#include "superposition/superposition_context.h"

#include <iostream>

namespace boost { namespace program_options { class variables_map; } }

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

using boost::lexical_cast;
using boost::ptr_vector;

const string cath_superpose_options::PROGRAM_NAME("cath-superpose");

/// TODOCUMENT
string cath_superpose_options::do_get_program_name() const {
	return PROGRAM_NAME;
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
string cath_superpose_options::do_update_error_or_help_string(const options_description &arg_visible_program_options ///< The full options_description of visible options
                                                              ) const {
	// If help was requested, then provide it
	if (the_misc_options_block.get_help()) {
		return the_misc_options_block.get_help_string(arg_visible_program_options, get_help_prefix_string(), get_help_suffix_string());
	}


	// If version information was requested, then provide it
	if (the_misc_options_block.get_version()) {
		return the_misc_options_block.get_version_string(get_program_name(), "This superposes protein structures.");
	}

	// Grab the objects from the options blocks
	const ptr_vector<alignment_acquirer> alignment_acquirers = the_alignment_input_options_block.get_alignment_acquirers();
	const ptr_vector<pdbs_acquirer>      pdb_acquirers       = the_pdb_input_options_block.get_pdbs_acquirers();
	const display_spec                   the_display_spec    = the_display_options_block.get_display_spec();
	const alignment_outputter_list       aln_outputters      = the_alignment_output_options_block.get_alignment_outputters( the_display_spec );
	const superposition_outputter_list   sup_outputters      = the_superposition_output_options_block.get_superposition_outputters( the_display_spec );

	// If there are no objects then no options were specified so just output the standard usage error string
	if (alignment_acquirers.empty() && pdb_acquirers.empty() && sup_outputters.empty()) {
		return get_standard_usage_error_string();
	}

	// Check that there is exactly one source of alignment
	if (alignment_acquirers.size() != 1) {
		return "Please specify one source of an alignment or superposition ("
		       + lexical_cast<string>(alignment_acquirers.size())
		       + " specified)\n"
		       + get_standard_usage_error_string();
	}

	// Check that there is exactly one source of PDBs
	if (pdb_acquirers.size() != 1) {
		return "Please specify one source of PDBs ("
		       + lexical_cast<string>(pdb_acquirers.size())
		       + " specified)\n"
		       + get_standard_usage_error_string();
	}

	const variables_map &vm = get_variables_map();

	const bool using_display_spec =    any_alignment_outputters_involve_display_spec    ( aln_outputters )
	                                || any_superposition_outputters_involve_display_spec( sup_outputters );
	if ( the_display_options_block.has_specified_options( vm ) && ! using_display_spec ) {
		return "Cannot specify display options because no display is being used in superposition/alignment outputters\n"
		       + get_standard_usage_error_string();
	}

	return "";
}

string cath_superpose_options::get_help_prefix_string() {
	ostringstream help_ss;
	help_ss << "Usage: " << PROGRAM_NAME << " alignment_source pdb_file_source [superposition_outputs]" << endl;
	help_ss << "Superpose structures based on a description of the superposition or an alignment"       << endl;
	help_ss << endl;
	help_ss << "Please specify:"                                                                        << endl;
	help_ss << " * one alignment"                                                                       << endl;
	help_ss << " * one method of reading PDB files (number to match the alignment)\n";
	return help_ss.str();
}

string cath_superpose_options::get_help_suffix_string() {
	ostringstream help_ss;
	help_ss << "Usage examples:"                                                                                                                << endl;
	help_ss << " * cath-superpose --ssap-aln-infile 1cukA1bvsA.list --pdb-infile $PDBDIR/1cukA --pdb-infile $PDBDIR/1bvsA --sup-to-pymol"       << endl;
	help_ss << "     (Superpose 1cukA and 1bvsA (in directory $PDBDIR) based on SSAP alignment file 1cukA1bvsA.list and then display in PyMOL)" << endl;
	help_ss << " * cat pdb1 end_file pdb2 end_file pdb3 | cath-superpose --pdbs-from-stdin --sup-to-stdout --res-name-align "                   << endl;
	help_ss << "     (Superpose the structures from stdin based on matching residue names and then write them to stdout [common Genome3D use case])\n";
	return help_ss.str();
}

/// TODOCUMENT
void cath_superpose_options::check_ok_to_use() const {
	if ( ! get_error_or_help_string().empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use invalid cath_superpose_options"));
	}
	if ( the_misc_options_block.get_help() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use cath_superpose_options that is set to do_nothing"));
	}
	const bool alignment_outputs_to_stdout     = the_alignment_output_options_block.outputs_to_stdout();
	const bool superposition_outputs_to_stdout = the_superposition_output_options_block.outputs_to_stdout();
	if ( alignment_outputs_to_stdout && superposition_outputs_to_stdout ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot output both alignment and superposition to stdout"));
	}
}

/// TODOCUMENT
unique_ptr<const pdbs_acquirer> cath_superpose_options::get_pdbs_acquirer() const {
	check_ok_to_use();
	const ptr_vector<pdbs_acquirer> pdb_acquirers = the_pdb_input_options_block.get_pdbs_acquirers();
	assert( pdb_acquirers.size() == 1 );
	return static_cast<unique_ptr<const pdbs_acquirer> >( pdb_acquirers.front().clone() );
}

/// \brief TODOCUMENT
selection_policy_acquirer cath_superpose_options::get_selection_policy_acquirer() const {
//	if ( the_alignment_input_options_block.get_residue_name_align() ) {
//		return selection_policy_acquirer( common_residue_select_all_policy(), common_atom_select_ca_policy() );
//	}
//	else {
//		return selection_policy_acquirer( common_residue_select_best_score_percent_policy(), common_atom_select_ca_policy() );
//	}
	return selection_policy_acquirer( common_residue_select_best_score_percent_policy(), common_atom_select_ca_policy() );
}

/// \brief TODOCUMENT
cath_superpose_options::cath_superpose_options() {
	super::add_options_block( the_misc_options_block                 );
	super::add_options_block( the_alignment_input_options_block      );
	super::add_options_block( the_ids_options_block                  );
	super::add_options_block( the_pdb_input_options_block            );
	super::add_options_block( the_alignment_output_options_block     );
	super::add_options_block( the_superposition_output_options_block );
	super::add_options_block( the_display_options_block               );
}

/// \brief TODOCUMENT
superposition_context cath_superpose_options::get_superposition_context(istream &arg_istream, ///< TODOCUMENT
                                                                        ostream &arg_stderr   ///< TODOCUMENT
                                                                        ) const {
	check_ok_to_use();
	const unique_ptr<const pdbs_acquirer> pdbs_acquirer_ptr = get_pdbs_acquirer();
	const pdb_list_str_vec_pair pdbs_and_names = pdbs_acquirer_ptr->get_pdbs_and_names( arg_istream, false );
	const pdb_list &raw_pdbs  = pdbs_and_names.first;
	const str_vec  &raw_names = pdbs_and_names.second;
	const str_vec   names     = ( the_ids_options_block.num_ids_specified() == raw_pdbs.size() ) ? get_all_ids( the_ids_options_block )
	                                                                                             : raw_names;
	const pdb_list  pdbs      = pdb_list_of_backbone_complete_subset_pdbs(raw_pdbs);

	if ( raw_pdbs.empty() ) {
		logger::log_and_exit(
			logger::return_code::NO_PDB_FILES_LOADED,
			"No valid PDBs were loaded"
		);
	}

	// For now, all superpositions come from alignment so throw if !acquires_alignment()
	const ptr_vector<alignment_acquirer> alignment_acquirers = the_alignment_input_options_block.get_alignment_acquirers();
	assert(alignment_acquirers.size() == 1); // This should already have been checked elsewhere
	const alignment_acquirer &the_alignment_acquirer = alignment_acquirers.front();
//	if (alignment_acquirers.size() > 1) {
//		BOOST_THROW_EXCEPTION(invalid_arg("Not yet implemented handlers to acquire a superposition without an alignment"));
//	}
//	if (!acquires_alignment()) {
//		BOOST_THROW_EXCEPTION(not_implemented_exception("Not yet implemented handlers to acquire a superposition without an alignment"));
//	}

	// Use the alignment_acquirer to get the alignment and corresponding spanning tree
	const pair<alignment, size_size_pair_vec> alignment_and_spanning_tree( the_alignment_acquirer.get_alignment_and_spanning_tree( pdbs ) );
	const alignment          &the_alignment = alignment_and_spanning_tree.first;
	const size_size_pair_vec &spanning_tree = alignment_and_spanning_tree.second;

	// \todo Remove this hacky code and fix it
	const path ssap_scores_file(the_alignment_input_options_block.get_ssap_scores_file());
	if ( ! ssap_scores_file.empty() ) {
		return superposition_context(
			raw_pdbs, ///< This should be raw_pdbs, not pdbs, so that superpositions include stripped residues (eg HETATM only residues). /// \todo Consider adding fast, simple test that ssap_scores_file superposition output includes HETATMs.
			names,
			hacky_multi_ssap_fuction(
				pdbs,
				names,
				spanning_tree,
				ssap_scores_file.parent_path(),
				get_selection_policy_acquirer(),
				arg_stderr
			),
			the_alignment
		);
	}

	// Construct an align_based_superposition_acquirer from the data and return the superposition it generates
	const align_based_superposition_acquirer aln_based_sup_acq(
		raw_pdbs,
		names,
		the_alignment,
		spanning_tree,
		get_selection_policy_acquirer()
	);
	return aln_based_sup_acq.get_superposition( arg_stderr );
}

/// \brief TODOCUMENT
alignment_outputter_list cath_superpose_options::get_alignment_outputters() const {
	const display_spec the_display_spec = the_display_options_block.get_display_spec();
	return the_alignment_output_options_block.get_alignment_outputters( the_display_spec );
}

/// \brief TODOCUMENT
superposition_outputter_list cath_superpose_options::get_superposition_outputters() const {
	const display_spec the_display_spec = the_display_options_block.get_display_spec();
	return the_superposition_output_options_block.get_superposition_outputters( the_display_spec );
}
