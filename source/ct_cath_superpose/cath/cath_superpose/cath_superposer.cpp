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

#include <optional>

#include <boost/ptr_container/ptr_vector.hpp>

#include <spdlog/spdlog.h>

#include "cath/acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "cath/acquirer/pdbs_acquirer/pdbs_acquirer.hpp"
#include "cath/acquirer/selection_policy_acquirer/selection_policy_acquirer.hpp"
#include "cath/acquirer/superposition_acquirer/align_based_superposition_acquirer.hpp"
#include "cath/alignment/alignment_context.hpp"
#include "cath/cath_superpose/options/cath_superpose_options.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/logger.hpp"
#include "cath/common/property_tree/read_from_json_file.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/outputter/alignment_outputter/alignment_outputter.hpp"
#include "cath/outputter/alignment_outputter/alignment_outputter_list.hpp"
#include "cath/outputter/superposition_outputter/superposition_outputter.hpp"
#include "cath/outputter/superposition_outputter/superposition_outputter_list.hpp"
#include "cath/superposition/superposition_context.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;
using namespace ::cath::sup;
using namespace ::std;

using ::std::nullopt;

/// \brief Perform a cath-superpose job as specified by the cath_superpose_options argument
///
/// The input and output stream parameters default to cin and cout respectively but are configurable,
/// primarily for testing purposes
void cath_superposer::superpose(const cath_superpose_options &prm_cath_superpose_options, ///< The details of the cath_superpose job to perform
                                istream                      &prm_istream,                ///< The istream from which any stdin-like input should be read
                                ostream                      &prm_stdout,                 ///< The ostream to which any stdout-like output should be written
                                ostream                      &prm_stderr                  ///< The ostream to which any stderr-like output should be written
                                ) {
	// If the options are invalid or specify to do_nothing, then just return
	const auto &error_or_help_string = prm_cath_superpose_options.get_error_or_help_string();
	if ( error_or_help_string ) {
		prm_stdout << *error_or_help_string;
		return;
	}

	// Get the superposition (and context) as specified by the cath_superpose_options
	const superposition_context the_sup_context = get_superposition_context( prm_cath_superpose_options, prm_istream, prm_stderr );

	// If there is an alignment, then for each of the alignment_outputters specified by the cath_superpose_options, output it
	const auto aln_outputters = prm_cath_superpose_options.get_alignment_outputters();
	if ( the_sup_context.has_alignment() ) {
		use_all_alignment_outputters(
			aln_outputters,
			make_restricted_alignment_context( the_sup_context ),
			prm_stdout,
			prm_stderr
		);
	}
	else {
		if ( ! aln_outputters.empty() ) {
			::spdlog::warn( "Ignoring alignment output options because there is no alignment (and the cake is a lie)" );
		}
	}

	// For each of the superposition_outputters specified by the cath_superpose_options, output the superposition
	const superposition_outputter_list outputters = prm_cath_superpose_options.get_superposition_outputters(
		aln_outputters.empty() ? default_supn_outputter::PYMOL : default_supn_outputter::NONE
	);
	use_all_superposition_outputters( outputters, the_sup_context, prm_stdout, prm_stderr );
}

/// \brief TODOCUMENT
///
/// \relates cath_superpose_options
///
/// \TODO Consider taking an ostream_ref_opt argument rather than ostream
///       (fix all errors, *then* provide default of ::std::nullopt)
superposition_context cath_superposer::get_superposition_context(const cath_superpose_options &prm_cath_sup_opts, ///< TODOCUMENT
                                                                 istream                      &prm_istream,       ///< TODOCUMENT
                                                                 ostream                      &prm_stderr         ///< TODOCUMENT
                                                                 ) {
	prm_cath_sup_opts.check_ok_to_use();

	// Grab the PDBs and their IDs
	const strucs_context  context  = get_pdbs_and_names( prm_cath_sup_opts, prm_istream, false );

	if ( context.get_pdbs().empty() ) {
		logger::log_and_exit(
			logger::return_code::NO_PDB_FILES_LOADED,
			"No valid PDBs were loaded"
		);
	}

	// If a JSON superposition file has been specified, return that
	const path_opt &json_sup_infile = prm_cath_sup_opts.get_json_sup_infile();
	if ( json_sup_infile ) {
		return set_pdbs_copy(
			read_from_json_file<superposition_context>( *json_sup_infile ),
			context.get_pdbs()
		);
	}

	// Get the alignment and corresponding spanning tree
	const auto aln_and_spn_tree = [&] {
		try {
			return get_alignment_and_spanning_tree(
				prm_cath_sup_opts,
				restrict_pdbs_copy( context ),
				get_align_refining( prm_cath_sup_opts )
				// ref( prm_stderr )
			);
		}
		catch (const std::exception &e) {
			logger::log_and_exit(
				logger::return_code::NO_PDB_FILES_LOADED,
				"Problem building alignment (and spanning tree) : "s + e.what()
			);
			exit( 0 ); //< Convince the compiler that this isn't a non-returning path out of the lambda
		}
	} ();
	const alignment &the_alignment = aln_and_spn_tree.first;
	const auto      &spanning_tree = aln_and_spn_tree.second;

	// Construct an align_based_superposition_acquirer from the data and return the superposition it generates
	const align_based_superposition_acquirer aln_based_sup_acq(
		the_alignment,
		spanning_tree,
		context,
		prm_cath_sup_opts.get_selection_policy_acquirer()
	);
	return aln_based_sup_acq.get_superposition( prm_stderr );
}
