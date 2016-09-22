/// \file
/// \brief The cath_score_align_options class definitions

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

#include "cath_score_align_options.h"

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

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath::align;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::lexical_cast;
using boost::ptr_vector;

const string cath_score_align_options::PROGRAM_NAME("cath-score-align");

/// TODOCUMENT
string cath_score_align_options::do_get_program_name() const {
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
string cath_score_align_options::do_update_error_or_help_string(const options_description &arg_visible_program_options ///< The full options_description of visible options
                                                                ) const {
	// If help was requested, then provide it
	if (the_misc_options_block.get_help()) {
		return the_misc_options_block.get_help_string(arg_visible_program_options, get_help_prefix_string(), get_help_suffix_string());
	}

	// If version information was requested, then provide it
	if (the_misc_options_block.get_version()) {
		return the_misc_options_block.get_version_string(get_program_name(), get_overview_string() );
	}

	// Grab the objects from the options blocks
	const ptr_vector<alignment_acquirer> alignment_acquirers = the_alignment_input_options_block.get_alignment_acquirers();
	const ptr_vector<pdbs_acquirer>      pdb_acquirers       = the_pdb_input_options_block.get_pdbs_acquirers();

	// If there are no objects then no options were specified so just output the standard usage error string
	if ( alignment_acquirers.empty() && pdb_acquirers.empty() ) {
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

	//const variables_map &vm = get_variables_map();

	return "";
}

string cath_score_align_options::get_help_prefix_string() {
	return "Usage: " + PROGRAM_NAME + " alignment_source protein_file_source [superposition_outputs]\n\n"
		+ get_overview_string() + R"(

Please specify:
 * one alignment
 * one method of reading proteins (number of proteins currently restricted to 2)";
}

string cath_score_align_options::get_help_suffix_string() {
	return "";
}

/// TODOCUMENT
void cath_score_align_options::check_ok_to_use() const {
	if (!get_error_or_help_string().empty()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use invalid cath_score_align_options"));
	}
	if (the_misc_options_block.get_help()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use cath_score_align_options that is set to do_nothing"));
	}
}


/// \brief TODOCUMENT
cath_score_align_options::cath_score_align_options() {
	super::add_options_block( the_misc_options_block            );
	super::add_options_block( the_alignment_input_options_block );
	super::add_options_block( the_pdb_input_options_block       );
}

/// TODOCUMENT
unique_ptr<const pdbs_acquirer> cath_score_align_options::get_pdbs_acquirer() const {
	check_ok_to_use();
	const ptr_vector<pdbs_acquirer> pdb_acquirers = the_pdb_input_options_block.get_pdbs_acquirers();
	assert( pdb_acquirers.size() == 1 );
	return static_cast<unique_ptr<const pdbs_acquirer> >( pdb_acquirers.front().clone() );
}

/// \brief TODOCUMENT
ptr_vector<alignment_acquirer> cath_score_align_options::get_alignment_acquirers() const {
	check_ok_to_use();
	return the_alignment_input_options_block.get_alignment_acquirers();
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string cath_score_align_options::get_overview_string() {
	return "Score an existing alignment using structural data";
}
