/// \file
/// \brief The cath_align_scorer class definitions

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

#include "cath_align_scorer.hpp"

#include "acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "acquirer/pdbs_acquirer/pdbs_acquirer.hpp"
#include "alignment/alignment.hpp"
#include "cath_score_align/options/cath_score_align_options.hpp"
#include "common/exception/not_implemented_exception.hpp"
#include "file/name_set/name_set_list.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_list_factory.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_value_list.hpp"
#include "score/aligned_pair_score_list/score_value_list_outputter/score_value_list_json_outputter.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::score;
using namespace std;

using boost::ptr_vector;

/// \brief Perform a cath-score-align job as specified by the cath_superpose_options argument
///
/// The input and output stream parameters default to cin and cout respectively but are configurable,
/// primarily for testing purposes
void cath_align_scorer::score(const cath_score_align_options &arg_cath_score_align_options, ///< The details of the cath-score-align job to perform
                              istream                         &arg_istream,                 ///< The istream from which any stdin-like input should be read
                              ostream                         &arg_stdout,                  ///< The ostream to which any stdout-like output should be written
                              ostream                         &/*arg_stderr*/               ///< The ostream to which any stderr-like output should be written
                              ) {
	// If the options are invalid or specify to do_nothing, then just return
	const auto &error_or_help_string = arg_cath_score_align_options.get_error_or_help_string();
	if ( error_or_help_string ) {
		arg_stdout << *error_or_help_string;
		return;
	}

	const auto           pdbs_acquirer_ptr = get_pdbs_acquirer( arg_cath_score_align_options );
	const auto           pdbs_and_names    = pdbs_acquirer_ptr->get_pdbs_and_names( arg_istream, true );
	const pdb_list      &pdbs              = pdbs_and_names.first;
	const name_set_list &names             = pdbs_and_names.second;
	const protein_list  proteins           = build_protein_list_of_pdb_list_and_names( pdbs, names );

	// An alignment is required but this should have been checked elsewhere
	const auto       alignment_acq_ptr  = get_alignment_acquirer( arg_cath_score_align_options );
	const auto       alignment_and_tree = alignment_acq_ptr->get_alignment_and_spanning_tree( pdbs );
	const alignment &the_alignment      = alignment_and_tree.first;
	// const auto      &spanning_tree      = alignment_and_tree.second;

	if ( proteins.size() != 2 || the_alignment.num_entries() != 2 ) {
		BOOST_THROW_EXCEPTION(not_implemented_exception("Currently only able to score alignments of more than two structures"));
	}

//	const alignment scored_scored_alignment = score_alignment_copy( residue_scorer(), scored_alignment, proteins );

	const protein &protein_a = proteins[ 0 ];
	const protein &protein_b = proteins[ 1 ];
	const aligned_pair_score_value_list the_scores_and_values = make_aligned_pair_score_value_list(

		make_default_aligned_pair_score_list(),
//		make_full_aligned_pair_score_list(),
//		make_seq_sim_flavours_aligned_pair_score_list(),
//		make_ssap_flavours_aligned_pair_score_list(),

		the_alignment,
		protein_a,
		protein_b
	);
	arg_stdout << score_value_list_json_outputter( the_scores_and_values ) << endl;
}
