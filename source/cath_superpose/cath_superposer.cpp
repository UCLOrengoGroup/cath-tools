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

#include "cath_superposer.hpp"

#include <boost/filesystem/path.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "acquirer/pdbs_acquirer/pdbs_acquirer.hpp"
#include "acquirer/selection_policy_acquirer/selection_policy_acquirer.hpp"
#include "acquirer/superposition_acquirer/align_based_superposition_acquirer.hpp"
#include "alignment/alignment_context.hpp"
#include "cath_superpose/options/cath_superpose_options.hpp"
#include "chopping/region/region.hpp"
#include "common/logger.hpp"
#include "common/property_tree/read_from_json_file.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "outputter/alignment_outputter/alignment_outputter.hpp"
#include "outputter/alignment_outputter/alignment_outputter_list.hpp"
#include "outputter/superposition_outputter/superposition_outputter.hpp"
#include "outputter/superposition_outputter/superposition_outputter_list.hpp"
#include "superposition/superposition_context.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
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
	const auto &error_or_help_string = arg_cath_superpose_options.get_error_or_help_string();
	if ( error_or_help_string ) {
		arg_stdout << *error_or_help_string << endl;
		return;
	}

	// Get the superposition (and context) as specified by the cath_superpose_options
	const superposition_context the_sup_context = get_superposition_context( arg_cath_superpose_options, arg_istream, arg_stderr );

	// If there is an alignment, then for each of the alignment_outputters specified by the cath_superpose_options, output it
	if ( the_sup_context.has_alignment() ) {
		use_all_alignment_outputters(
			arg_cath_superpose_options.get_alignment_outputters(),
			make_alignment_context( the_sup_context ),
			arg_stdout,
			arg_stderr
		);
	}

	// For each of the superposition_outputters specified by the cath_superpose_options, output the superposition
	const superposition_outputter_list outputters = arg_cath_superpose_options.get_superposition_outputters();
	use_all_superposition_outputters( outputters, the_sup_context, arg_stdout, arg_stderr );
}

/// \brief TODOCUMENT
///
/// \relates cath_superpose_options
///
/// \TODO Consider taking an ostream_ref_opt argument rather than ostream
///       (fix all errors, *then* provide default of boost::none)
superposition_context cath_superposer::get_superposition_context(const cath_superpose_options &arg_cath_sup_opts, ///< TODOCUMENT
                                                                 istream                      &arg_istream,       ///< TODOCUMENT
                                                                 ostream                      &arg_stderr         ///< TODOCUMENT
                                                                 ) {
	arg_cath_sup_opts.check_ok_to_use();

	// Grab the PDBs and their IDs
	const strucs_context  context  = get_pdbs_and_names( arg_cath_sup_opts, arg_istream, false );
	const pdb_list       &raw_pdbs = context.get_pdbs();

	if ( raw_pdbs.empty() ) {
		logger::log_and_exit(
			logger::return_code::NO_PDB_FILES_LOADED,
			"No valid PDBs were loaded"
		);
	}

	// If a JSON superposition file has been specified, return that
	const path_opt &json_sup_infile = arg_cath_sup_opts.get_json_sup_infile();
	if ( json_sup_infile ) {
		return set_pdbs_copy(
			read_from_json_file<superposition_context>( *json_sup_infile ),
			raw_pdbs
		);
	}

	// Get the alignment and corresponding spanning tree
	const pdb_list backbone_complete_subset_pdbs = pdb_list_of_backbone_complete_region_limited_subset_pdbs(
		raw_pdbs,
		context.get_regions(),
		ref( arg_stderr )
	);
	const auto       aln_and_spn_tree = [&] () {
		try {
			return get_alignment_and_spanning_tree( arg_cath_sup_opts, backbone_complete_subset_pdbs );
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
	const path ssap_scores_file = arg_cath_sup_opts.get_alignment_input_spec().get_ssap_scores_file();
	if ( ! ssap_scores_file.empty() ) {
		return {
			hacky_multi_ssap_fuction(
				backbone_complete_subset_pdbs,
				get_multi_ssap_alignment_file_names( context.get_name_sets() ),
				spanning_tree,
				ssap_scores_file.parent_path(),
				arg_cath_sup_opts.get_selection_policy_acquirer(),
				arg_stderr
			),
			context, ///< Importantly, this contains the raw_pdbs, not backbone_complete_subset_pdbs, so that superpositions include stripped residues (eg HETATM only residues). /// \todo Consider adding fast, simple test that ssap_scores_file superposition output includes HETATMs.
			the_alignment
		};
	}

	// Construct an align_based_superposition_acquirer from the data and return the superposition it generates
	const align_based_superposition_acquirer aln_based_sup_acq(
		the_alignment,
		spanning_tree,
		context,
		arg_cath_sup_opts.get_selection_policy_acquirer()
	);
	return aln_based_sup_acq.get_superposition( arg_stderr );
}
