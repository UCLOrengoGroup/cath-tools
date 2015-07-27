/// \file
/// \brief The cath_refine_align_options class definitions

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

#include "cath_refine_align_options.h"

#include <boost/program_options.hpp>
#include <boost/shared_array.hpp>

#include "alignment/alignment.h"
#include "alignment/common_atom_selection_policy/common_atom_select_ca_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.h"
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
#include "options/outputter/alignment_outputter/alignment_outputter.h"
#include "options/outputter/alignment_outputter/alignment_outputter_list.h"
#include "options/outputter/superposition_outputter/superposition_outputter.h"
#include "options/outputter/superposition_outputter/superposition_outputter_list.h"
#include "superposition/superposition_context.h"

#include <iostream>

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath::align;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::lexical_cast;
using boost::ptr_vector;

const string cath_refine_align_options::PROGRAM_NAME("cath-refine-align");

/// TODOCUMENT
string cath_refine_align_options::do_get_program_name() const {
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
string cath_refine_align_options::do_update_error_or_help_string(const options_description &arg_visible_program_options ///< The full options_description of visible options
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

string cath_refine_align_options::get_help_prefix_string() {
	ostringstream help_ss;
	help_ss << "Usage: " << PROGRAM_NAME << " alignment_source protein_file_source [superposition_outputs]" << endl;
	help_ss << "Refine an existing alignment by iteratively attempting to optimise SSAP score"          << endl;
	help_ss << endl;
	help_ss << "Please specify:"                                                                        << endl;
	help_ss << " * one alignment"                                                                       << endl;
	help_ss << " * one method of reading proteins (number of proteins currently restricted to 2)\n";
	return help_ss.str();
}

string cath_refine_align_options::get_help_suffix_string() {
	ostringstream help_ss;
	help_ss << "Usage examples:"                                                                                                                << endl;
	help_ss << " * cath-superpose --ssap-aln-infile 1cukA1bvsA.list --pdb-infile $PDBDIR/1cukA --pdb-infile $PDBDIR/1bvsA --sup-to-pymol"       << endl;
	help_ss << "     (Superpose 1cukA and 1bvsA (in directory $PDBDIR) based on SSAP alignment file 1cukA1bvsA.list and then display in PyMOL)" << endl;
	help_ss << " * cat pdb1 end_file pdb2 end_file pdb3 | cath-superpose --pdbs-from-stdin --sup-to-stdout --res-name-align "                   << endl;
	help_ss << "     (Superpose the structures from stdin based on matching residue names and then write them to stdout [common Genome3D use case])\n";
	return help_ss.str();
}

/// TODOCUMENT
void cath_refine_align_options::check_ok_to_use() const {
	if (!get_error_or_help_string().empty()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use invalid cath_refine_align_options"));
	}
	if (the_misc_options_block.get_help()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use cath_refine_align_options that is set to do_nothing"));
	}
	const bool alignment_outputs_to_stdout     = the_alignment_output_options_block.outputs_to_stdout();
	const bool superposition_outputs_to_stdout = the_superposition_output_options_block.outputs_to_stdout();
	if ( alignment_outputs_to_stdout && superposition_outputs_to_stdout ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot output both alignment and superposition to stdout"));
	}
}


/// \brief TODOCUMENT
cath_refine_align_options::cath_refine_align_options() {
	super::add_options_block(the_misc_options_block);
	super::add_options_block(the_alignment_input_options_block);
	super::add_options_block(the_pdb_input_options_block);
	super::add_options_block(the_alignment_output_options_block);
	super::add_options_block(the_superposition_output_options_block);
	super::add_options_block(the_display_options_block);
}

/// TODOCUMENT
unique_ptr<const pdbs_acquirer> cath_refine_align_options::get_pdbs_acquirer() const {
	check_ok_to_use();
	const ptr_vector<pdbs_acquirer> pdb_acquirers = the_pdb_input_options_block.get_pdbs_acquirers();
	assert( pdb_acquirers.size() == 1 );
	return static_cast<unique_ptr<const pdbs_acquirer> >( pdb_acquirers.front().clone() );
}

/// \brief TODOCUMENT
ptr_vector<alignment_acquirer> cath_refine_align_options::get_alignment_acquirers() const {
	check_ok_to_use();
	return the_alignment_input_options_block.get_alignment_acquirers();
}

/// \brief TODOCUMENT
selection_policy_acquirer cath_refine_align_options::get_selection_policy_acquirer() const {
	check_ok_to_use();
	if ( the_alignment_input_options_block.get_residue_name_align() ) {
		return selection_policy_acquirer( common_residue_select_all_policy(), common_atom_select_ca_policy() );
	}
	else {
		return selection_policy_acquirer( common_residue_select_best_score_percent_policy(), common_atom_select_ca_policy() );
	}
}

/// \brief TODOCUMENT
alignment_outputter_list cath_refine_align_options::get_alignment_outputters() const {
	check_ok_to_use();
	const display_spec the_display_spec = the_display_options_block.get_display_spec();
	return the_alignment_output_options_block.get_alignment_outputters( the_display_spec );
}

/// \brief TODOCUMENT
superposition_outputter_list cath_refine_align_options::get_superposition_outputters() const {
	check_ok_to_use();
	const display_spec the_display_spec = the_display_options_block.get_display_spec();
	return the_superposition_output_options_block.get_superposition_outputters( the_display_spec );
}
