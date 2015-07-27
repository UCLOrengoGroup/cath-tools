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

#include "cath_align_scorer.h"

#include "alignment/alignment.h"
//#include "alignment/alignment_context.h"
//#include "alignment/gap/gap_penalty.h"
//#include "alignment/residue_score/residue_scorer.h"
#include "exception/not_implemented_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "file/pdb/pdb_list.h"
#include "options/acquirer/alignment_acquirer/alignment_acquirer.h"
#include "options/acquirer/pdbs_acquirer/pdbs_acquirer.h"
#include "options/executable/cath_score_align_options/cath_score_align_options.h"
#include "score/aligned_pair_score_list/aligned_pair_score_list_factory.h"
#include "score/aligned_pair_score_list/aligned_pair_score_value_list.h"
#include "score/aligned_pair_score_list/score_value_list_outputter/score_value_list_json_outputter.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_list.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

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
	const string error_or_help_string = arg_cath_score_align_options.get_error_or_help_string();
	if (!error_or_help_string.empty()) {
		arg_stdout << error_or_help_string << endl;
		return;
	}

	const unique_ptr<const pdbs_acquirer> pdbs_acquirer_ptr = arg_cath_score_align_options.get_pdbs_acquirer();
	const pdb_list_str_vec_pair           pdbs_and_names    = pdbs_acquirer_ptr->get_pdbs_and_names( arg_istream, true );
	const pdb_list                       &pdbs              = pdbs_and_names.first;
	const str_vec                        &names             = pdbs_and_names.second;
	const protein_list                    proteins          = build_protein_list_of_pdb_list_and_names( pdbs, names );

	// An alignment is required but this should have been checked elsewhere
	const ptr_vector<alignment_acquirer> alignment_acquirers = arg_cath_score_align_options.get_alignment_acquirers();
	assert(alignment_acquirers.size() == 1); // This should already have been checked elsewhere
	const alignment_acquirer &the_alignment_acquirer = alignment_acquirers.front();

	const pair<alignment, size_size_pair_vec>  alignment_and_tree = the_alignment_acquirer.get_alignment_and_spanning_tree(pdbs);
	const alignment                           &the_alignment      = alignment_and_tree.first;
	// const size_size_pair_vec                  &spanning_tree      = alignment_and_tree.second;

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
