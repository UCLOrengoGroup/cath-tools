/// \file
/// \brief The cath_superposer class definitions

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

#include "cath_superposer.h"

#include <boost/filesystem/path.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "acquirer/alignment_acquirer/alignment_acquirer.h"
#include "acquirer/pdbs_acquirer/pdbs_acquirer.h"
#include "acquirer/selection_policy_acquirer/selection_policy_acquirer.h"
#include "acquirer/superposition_acquirer/align_based_superposition_acquirer.h"
#include "alignment/alignment_context.h"
#include "cath_superpose/options/cath_superpose_options.h"
#include "common/logger.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "outputter/alignment_outputter/alignment_outputter.h"
#include "outputter/alignment_outputter/alignment_outputter_list.h"
#include "outputter/superposition_outputter/superposition_outputter.h"
#include "outputter/superposition_outputter/superposition_outputter_list.h"
#include "superposition/superposition_context.h"

using namespace cath;
using namespace cath::align;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

using boost::filesystem::path;
using boost::ptr_vector;

/// \brief Perform a cath-superpose job as specified by the cath_superpose_options argument
///
/// The input and output stream parameters default to cin and cout respectively but are configurable,
/// primarily for testing purposes
void cath_superposer::superpose(const cath_superpose_options &arg_cath_superpose_options, ///< The details of the cath_superpose job to perform
                                istream                      &arg_istream,                ///< The istream from which any stdin-like input should be read
                                ostream                      &arg_stdout,                 ///< The ostream to which any stdout-like output should be written
                                ostream                      &arg_stderr                  ///< The ostream to which any stderr-like output should be written
                                ) {
	// If the options are invalid or specify to do_nothing, then just return
	const auto error_or_help_string = arg_cath_superpose_options.get_error_or_help_string();
	if ( error_or_help_string ) {
		arg_stderr << *error_or_help_string << endl;
		return;
	}

	// Get the superposition (and context) as specified by the cath_superpose_options
	const superposition_context the_sup_context = get_superposition_context( arg_cath_superpose_options, arg_istream, arg_stderr );

	// If there is an alignment, then for each of the alignment_outputters specified by the cath_superpose_options, output it
	if ( the_sup_context.has_alignment() ) {
		const pdb_list                 &pdbs           = the_sup_context.get_pdbs_cref();
		const str_vec                  &names          = the_sup_context.get_names_cref();
		const alignment                &the_aln        = the_sup_context.get_alignment_cref();
		const alignment_outputter_list  aln_outputters = arg_cath_superpose_options.get_alignment_outputters();
		use_all_alignment_outputters( aln_outputters, alignment_context( pdbs, names, the_aln ), arg_stdout, arg_stderr );
	}

	// For each of the superposition_outputters specified by the cath_superpose_options, output the superposition
	const superposition_outputter_list outputters = arg_cath_superpose_options.get_superposition_outputters();
	use_all_superposition_outputters( outputters, the_sup_context, arg_stdout, arg_stderr );
}

/// \brief TODOCUMENT
///
/// \relates cath_superpose_options
superposition_context cath_superposer::get_superposition_context(const cath_superpose_options &arg_cath_superpose_options, ///< TODOCUMENT
                                                                 istream                      &arg_istream,                ///< TODOCUMENT
                                                                 ostream                      &arg_stderr                  ///< TODOCUMENT
                                                                 ) {
	arg_cath_superpose_options.check_ok_to_use();

	// Grab the PDBs and their IDs
	const auto      pdbs_and_names = get_pdbs_and_names( arg_cath_superpose_options, arg_istream, false );
	const pdb_list &raw_pdbs       = pdbs_and_names.first;
	const str_vec  &raw_names      = pdbs_and_names.second;
	const str_vec  &ids_block_ids  = arg_cath_superpose_options.get_ids();
	const str_vec  &names          = ( ids_block_ids.size() == raw_pdbs.size() ) ? ids_block_ids
	                                                                             : raw_names;
	const pdb_list  pdbs           = pdb_list_of_backbone_complete_subset_pdbs( raw_pdbs );

	if ( raw_pdbs.empty() ) {
		logger::log_and_exit(
			logger::return_code::NO_PDB_FILES_LOADED,
			"No valid PDBs were loaded"
		);
	}

	// If a JSON superposition file has ben specified, return that
	const path_opt json_sup_infile = arg_cath_superpose_options.get_json_sup_infile();
	if ( json_sup_infile ) {
		return set_pdbs_copy(
			read_superposition_context_from_json_file( *json_sup_infile ),
			pdbs
		);
	}

	// Get the alignment and corresponding spanning tree
	const auto       aln_and_spn_tree = [&] () {
		try {
			return get_alignment_and_spanning_tree( arg_cath_superpose_options, pdbs );
		}
		catch (const std::exception &e) {
			logger::log_and_exit(
				logger::return_code::NO_PDB_FILES_LOADED,
				"Problem building alignment (and spanning tree) : "s + e.what()
			);
			exit( 0 ); //< Convince the compiler that this isn't a non-returning path out of the lambda
		}
	}();
	const alignment &the_alignment    = aln_and_spn_tree.first;
	const auto      &spanning_tree    = aln_and_spn_tree.second;

	// \todo Remove this hacky code and fix it
	const path ssap_scores_file = arg_cath_superpose_options.get_alignment_input_spec().get_ssap_scores_file();
	if ( ! ssap_scores_file.empty() ) {
		return superposition_context(
			raw_pdbs, ///< This should be raw_pdbs, not pdbs, so that superpositions include stripped residues (eg HETATM only residues). /// \todo Consider adding fast, simple test that ssap_scores_file superposition output includes HETATMs.
			names,
			hacky_multi_ssap_fuction(
				pdbs,
				names,
				spanning_tree,
				ssap_scores_file.parent_path(),
				arg_cath_superpose_options.get_selection_policy_acquirer(),
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
		arg_cath_superpose_options.get_selection_policy_acquirer()
	);
	return aln_based_sup_acq.get_superposition( arg_stderr );
}
