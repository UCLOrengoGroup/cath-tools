/// \file
/// \brief The cath_align_refiner class definitions

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

#include "cath_align_refiner.hpp"

#include "acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "acquirer/pdbs_acquirer/pdbs_acquirer.hpp"
#include "acquirer/superposition_acquirer/align_based_superposition_acquirer.hpp"
#include "alignment/alignment.hpp"
#include "alignment/alignment_context.hpp"
#include "alignment/gap/gap_penalty.hpp"
#include "alignment/refiner/alignment_refiner.hpp"
#include "alignment/residue_score/residue_scorer.hpp"
#include "cath_refine_align/options/cath_refine_align_options.hpp"
#include "chopping/region/region.hpp"
#include "common/exception/not_implemented_exception.hpp"
#include "outputter/alignment_outputter/alignment_outputter.hpp"
#include "outputter/alignment_outputter/alignment_outputter_list.hpp"
#include "outputter/superposition_outputter/superposition_outputter.hpp"
#include "outputter/superposition_outputter/superposition_outputter_list.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_list_factory.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_value_list.hpp"
#include "score/aligned_pair_score_list/score_value_list_outputter/score_value_list_json_outputter.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "superposition/superposition_context.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::align::gap;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::score;
using namespace std;

/// \brief Perform a cath-refine-align job as specified by the cath_superpose_options argument
///
/// The input and output stream parameters default to cin and cout respectively but are configurable,
/// primarily for testing purposes
void cath_align_refiner::refine(const cath_refine_align_options &arg_cath_refine_align_options, ///< The details of the cath-refine-align job to perform
                                istream                         &arg_istream,                   ///< The istream from which any stdin-like input should be read
                                ostream                         &arg_stdout,                   ///< The ostream to which any stdout-like output should be written
                                ostream                         &arg_stderr                     ///< The ostream to which any stderr-like output should be written
                                ) {
	// If the options are invalid or specify to do_nothing, then just return
	const auto &error_or_help_string = arg_cath_refine_align_options.get_error_or_help_string();
	if ( error_or_help_string ) {
		arg_stdout << *error_or_help_string;
		return;
	}

	// Grab the PDBs and their IDs
	const strucs_context  context      = get_pdbs_and_names( arg_cath_refine_align_options, arg_istream, true );

	// TODO: Populate this name summarising the overall group so it can be used, eg in superposition scripts
	const string          overall_name{};

	const auto backbone_complete_strucs_context = strucs_context_of_backbone_complete_region_limited_subset_pdbs(
		context,
		ref( arg_stderr )
	);

	// An alignment is required but this should have been checked elsewhere
	const protein_list  proteins           = build_protein_list( backbone_complete_strucs_context );
	const auto          aln_acq_ptr        = get_alignment_acquirer( arg_cath_refine_align_options );
	const auto          alignment_and_tree = aln_acq_ptr->get_alignment_and_spanning_tree( backbone_complete_strucs_context );
	const alignment    &the_alignment      = alignment_and_tree.first;
	const auto         &spanning_tree      = alignment_and_tree.second;

//	if ( proteins.size() != 2 || the_alignment.num_entries() != 2 ) {
//		BOOST_THROW_EXCEPTION(not_implemented_exception("Currently only able to score alignments of more than two structures"));
//	}
//	const protein &protein_a = proteins[0];
//	const protein &protein_b = proteins[1];
//	const aligned_pair_score_value_list the_scores_and_values = make_aligned_pair_score_value_list(
//		make_full_aligned_pair_score_list(),
//		the_alignment,
//		protein_a,
//		protein_b
//	);
////	cout << "Scores are as follows : " << endl;
//	cout << the_scores_and_values << flush;
//
//	return;

	const alignment refined_alignment        = alignment_refiner().iterate( the_alignment, proteins, gap_penalty( 50, 0 ) );
	const alignment scored_refined_alignment = score_alignment_copy( residue_scorer(), refined_alignment, proteins );

//	const protein &protein_a = proteins[0];
//	const protein &protein_b = proteins[1];
//	const aligned_pair_score_value_list the_scores_and_values = make_aligned_pair_score_value_list(
//		make_full_aligned_pair_score_list(),
//		scored_refined_alignment,
//		protein_a,
//		protein_b
//	);
//	cout << "Scores are as follows : " << endl;
//	cout << score_value_list_json_outputter( the_scores_and_values ) << endl;
////	cout << the_scores_and_values << flush;
////	write_alignment_as_fasta_alignment( cout, refined_alignment, proteins );

	// Construct an align_based_superposition_acquirer from the data and return the superposition it generates
	const align_based_superposition_acquirer aln_based_sup_acq(
		scored_refined_alignment,
		spanning_tree,
		context,
		get_selection_policy_acquirer( arg_cath_refine_align_options )
	);

	// For each of the alignment_outputters specified by the cath_superpose_options, output the alignment
	const alignment_outputter_list aln_outputters = arg_cath_refine_align_options.get_alignment_outputters();
	use_all_alignment_outputters( aln_outputters, alignment_context{ scored_refined_alignment, context }, arg_stdout, arg_stderr );

	// For each of the superposition_outputters specified by the cath_superpose_options, output the superposition
	const superposition_outputter_list sup_outputters = arg_cath_refine_align_options.get_superposition_outputters(
		aln_outputters.empty() ? default_supn_outputter::PYMOL : default_supn_outputter::NONE
	);
	use_all_superposition_outputters( sup_outputters, aln_based_sup_acq.get_superposition( arg_stderr ), arg_stdout, arg_stderr, overall_name );
}
