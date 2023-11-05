/// \file
/// \brief The alignment_refiner class definitions

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

#include "alignment_refiner.hpp"

#include <fstream>

#include <boost/numeric/conversion/cast.hpp>

#include <spdlog/spdlog.h>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/dyn_prog_align/detail/matrix_plotter/gnuplot_matrix_plotter.hpp"
#include "cath/alignment/dyn_prog_align/detail/matrix_plotter/matrix_plot.hpp"
#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/new_matrix_dyn_prog_score_source.hpp"
#include "cath/alignment/dyn_prog_align/ssap_code_dyn_prog_aligner.hpp" // ***** TEMPORARY *****
#include "cath/alignment/dyn_prog_align/std_dyn_prog_aligner.hpp"
#include "cath/alignment/io/alignment_io.hpp"
#include "cath/alignment/pair_alignment.hpp"
#include "cath/alignment/refiner/detail/alignment_split.hpp"
#include "cath/alignment/refiner/detail/alignment_split_list.hpp"
#include "cath/alignment/refiner/detail/alignment_split_mapping.hpp"
#include "cath/alignment/residue_score/residue_scorer.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/ssap/ssap.hpp"
#include "cath/structure/entry_querier/residue_querier.hpp" // ***** TEMPORARY *****
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_list.hpp"
#include "cath/structure/view_cache/view_cache.hpp"
#include "cath/structure/view_cache/view_cache_list.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::align::detail;
using namespace ::cath::align::gap;
using namespace ::cath::common;
using namespace ::cath::index;
using namespace ::std;

/// \brief TODOCUMENT
bool_aln_pair alignment_refiner::iterate_step(const alignment       &prm_alignment,       ///< TODOCUMENT
                                              const protein_list    &prm_proteins,        ///< TODOCUMENT
                                              const view_cache_list &prm_view_cache_list, ///< TODOCUMENT
                                              const gap_penalty     &prm_gap_penalty      ///< TODOCUMENT
                                              ) {
	::spdlog::info( "Will search for sensible ways to split alignment with {} entries", prm_alignment.num_entries() );

	const size_t num_entries = prm_alignment.num_entries();
	if ( prm_proteins.size() != num_entries ) {
		BOOST_THROW_EXCEPTION(not_implemented_exception("Mismatch between number of entries in alignment and in protein list"));
	}

	return iterate_step_for_alignment_split_list(
		prm_alignment,
		prm_proteins,
		prm_view_cache_list,
		prm_gap_penalty,
		get_standard_alignment_splits( prm_alignment )
	);
}

/// \brief TODOCUMENT
bool_aln_pair alignment_refiner::iterate_step_for_alignment_split_list(const alignment            &prm_alignment,           ///< TODOCUMENT
                                                                       const protein_list         &prm_proteins,            ///< TODOCUMENT
                                                                       const view_cache_list      &prm_view_cache_list,     ///< TODOCUMENT
                                                                       const gap_penalty          &prm_gap_penalty,         ///< TODOCUMENT
                                                                       const alignment_split_list &prm_alignment_split_list ///< TODOCUMENT
                                                                       ) {
	alignment iter_aln( prm_alignment );
	bool inserted_residues = false;
	for (const alignment_split &the_split : prm_alignment_split_list) {
		const bool_aln_pair inserted_res_and_aln = iterate_step_for_alignment_split(
			iter_aln,
			prm_proteins,
			prm_view_cache_list,
			prm_gap_penalty,
			the_split
		);
		inserted_residues = ( inserted_res_and_aln.first || inserted_residues );
		iter_aln          =   inserted_res_and_aln.second;
	}
	return make_pair( inserted_residues, iter_aln );
}

/// \brief TODOCUMENT
bool_aln_pair alignment_refiner::iterate_step_for_alignment_split(const alignment       &prm_alignment,       ///< TODOCUMENT
                                                                  const protein_list    &prm_proteins,        ///< TODOCUMENT
                                                                  const view_cache_list &prm_view_cache_list, ///< TODOCUMENT
                                                                  const gap_penalty     &prm_gap_penalty,     ///< TODOCUMENT
                                                                  const alignment_split &prm_alignment_split  ///< TODOCUMENT
                                                                  ) {
	const size_vec correct_lengths = get_protein_lengths( prm_proteins );

	::spdlog::info( "Iterating alignment with {} entries", prm_alignment.num_entries() );
	// cerr << "Iterating alignment with " << prm_alignment.num_entries() << " entries using split :";
	// for (const size_t &split_member : prm_alignment_split) {
	// 	cerr << " " << split_member;
	// }
	// cerr << endl;

	if ( prm_alignment_split.get_num_entries() != prm_alignment.num_entries() ) {
		BOOST_THROW_EXCEPTION(not_implemented_exception("Number of entries in alignment split doesn't match number in alignment"));
	}

	// bob
	const alignment_split_mapping mapping_a = make_alignment_split_mapping( prm_alignment, prm_alignment_split, alignment_split_half::FIRST,  correct_lengths );
	const alignment_split_mapping mapping_b = make_alignment_split_mapping( prm_alignment, prm_alignment_split, alignment_split_half::SECOND, correct_lengths );
	const bool inserted_residues = ( mapping_a.inserted_entries() || mapping_b.inserted_entries() );

	// Grab the lengths
	const size_t full_length_a     = mapping_a.length();
	const size_t full_length_b     = mapping_b.length();
	const size_t full_window_width = get_window_width_for_full_matrix( full_length_a, full_length_b );

	from_alignment_scores.assign( full_length_a, float_score_vec( full_length_b, 0.0 ) );
	to_alignment_scores.assign  ( full_length_a, float_score_vec( full_length_b, 0.0 ) );

//	cerr << "number of entries in alignment is " << prm_alignment.num_entries() << endl;
//	cerr << "number of entries in half a of split is " << mapping_a.num_entries() << endl;
//	cerr << "number of entries in half b of split is " << mapping_b.num_entries() << endl;

	const alignment::size_type alignment_length = prm_alignment.length();
	for (const size_t &aln_ctr : indices( alignment_length ) ) {
		const size_opt mapping_index_a = mapping_a.index_of_orig_aln_index( aln_ctr );
		const size_opt mapping_index_b = mapping_b.index_of_orig_aln_index( aln_ctr );
		if ( ! mapping_index_a || ! mapping_index_b ) {
			continue;
		}
//		cerr << "At alignment counter " << aln_ctr << endl;

		const size_vec present_orig_aln_entries_a = present_orig_aln_entries_of_index( mapping_a, *mapping_index_a );
		const size_vec present_orig_aln_entries_b = present_orig_aln_entries_of_index( mapping_b, *mapping_index_b );

		for (const size_t &present_orig_aln_entry_a : present_orig_aln_entries_a) {
			for (const size_t &present_orig_aln_entry_b : present_orig_aln_entries_b) {
				const size_t        present_entry_a = * mapping_a.entry_of_orig_aln_entry( present_orig_aln_entry_a );
				const size_t        present_entry_b = * mapping_b.entry_of_orig_aln_entry( present_orig_aln_entry_b );
				const aln_posn_type a_position      = get_position_of_entry_of_index( mapping_a, present_entry_a, *mapping_index_a );
				const aln_posn_type b_position      = get_position_of_entry_of_index( mapping_b, present_entry_b, *mapping_index_b );

	//			cerr << "Rescoring based on pair " << a_position << ", " << b_position << endl;
				const protein &protein_a = prm_proteins[ present_orig_aln_entry_a ];
				const protein &protein_b = prm_proteins[ present_orig_aln_entry_b ];
				const size_t   length_a  = protein_a.get_length();
				const size_t   length_b  = protein_b.get_length();
//				const residue &residue_a = protein_a.get_residue_ref_of_index( a_position );
//				const residue &residue_b = protein_b.get_residue_ref_of_index( b_position );

				for (const size_t &res_ctr_a : indices( length_a ) ) {
					for (const size_t &res_ctr_b : indices( length_b ) ) {
						if ( res_ctr_a != a_position && res_ctr_b != b_position ) {
							const size_t   other_mapping_index_a = mapping_a.index_of_protein_index( present_entry_a, res_ctr_a );
							const size_t   other_mapping_index_b = mapping_b.index_of_protein_index( present_entry_b, res_ctr_b );
//							const residue &other_residue_a       = protein_a.get_residue_ref_of_index( res_ctr_a );
//							const residue &other_residue_b       = protein_b.get_residue_ref_of_index( res_ctr_b );

							from_alignment_scores[ other_mapping_index_a ][ other_mapping_index_b ] += get_residue_context(
								prm_view_cache_list,
								present_orig_aln_entry_a,
								present_orig_aln_entry_b,
								a_position,
								b_position,
								res_ctr_a,
								res_ctr_b
							);
							to_alignment_scores[ other_mapping_index_a ][ other_mapping_index_b ] += get_residue_context(
								prm_view_cache_list,
								present_orig_aln_entry_a,
								present_orig_aln_entry_b,
								res_ctr_a,
								res_ctr_b,
								a_position,
								b_position
							);
						}
					}
				}
			}
		}
	}

	float_score_vec_vec avg_scores( full_length_a, float_score_vec( full_length_b, 0 ) );
	for (const size_t &ctr_a : indices( full_length_a ) ) {
		for (const size_t &ctr_b : indices( full_length_b ) ) {
			avg_scores[ctr_a][ctr_b] = (from_alignment_scores[ctr_a][ctr_b] + to_alignment_scores[ctr_a][ctr_b]) / 2.0;
		}
	}

//	matrix_plot<gnuplot_matrix_plotter>( path("matrix_plot_from"), new_matrix_dyn_prog_score_source( from_alignment_scores, length_a, length_b) );
//	matrix_plot<gnuplot_matrix_plotter>( path("matrix_plot___to"), new_matrix_dyn_prog_score_source( to_alignment_scores,   length_a, length_b) );
//	matrix_plot<gnuplot_matrix_plotter>( path("matrix_plot__avg"), new_matrix_dyn_prog_score_source( avg_scores,            length_a, length_b) );

	const new_matrix_dyn_prog_score_source scorer( avg_scores, full_length_a, full_length_b );

	const score_alignment_pair score_and_alignment = std_dyn_prog_aligner().align( scorer, prm_gap_penalty, full_window_width );

	alignment new_alignment = set_empty_scores_copy(
		build_alignment(
			score_and_alignment.second,
			mapping_a,
			mapping_b
		)
	);
	return make_pair( inserted_residues, new_alignment );
}

/// \brief TODOCUMENT
alignment alignment_refiner::iterate(const alignment    &prm_alignment,  ///< TODOCUMENT
                                     const protein_list &prm_proteins,   ///< TODOCUMENT
                                     const gap_penalty  &prm_gap_penalty ///< TODOCUMENT
                                     ) {
	return iterate(
		prm_alignment,
		prm_proteins,
		view_cache_list( prm_proteins ),
		prm_gap_penalty
	);
}

/// \brief TODOCUMENT
alignment alignment_refiner::iterate(const alignment       &prm_alignment,       ///< TODOCUMENT
                                     const protein_list    &prm_proteins,        ///< TODOCUMENT
                                     const view_cache_list &prm_view_cache_list, ///< TODOCUMENT
                                     const gap_penalty     &prm_gap_penalty      ///< TODOCUMENT
                                     ) {
//	if (prm_proteins.size() != 2 || prm_alignment.num_entries() != 2) {
//		BOOST_THROW_EXCEPTION(not_implemented_exception("Currently only able to iterate alignments of more than two structures"));
//	}

	/// \todo Move to using scores to prevent loops

	/// \todo Ensure that if using loops, a step that fills in alignment holes is always accepted

	// size_t iter_ctr = 0;
	alignment prev_alignment( alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT );
	alignment curr_alignment( prm_alignment );
	bool inserted_residues = true;
	while ( inserted_residues || curr_alignment != prev_alignment ) {
		const bool_aln_pair ins_res_and_next_aln = iterate_step( curr_alignment, prm_proteins, prm_view_cache_list, prm_gap_penalty );
		inserted_residues        = ins_res_and_next_aln.first;
		alignment next_alignment = ins_res_and_next_aln.second;
		const bool next_matches_prev= ( next_alignment == prev_alignment );
		swap( prev_alignment, curr_alignment );
		swap( curr_alignment, next_alignment );

		// if (iter_ctr > 0) {
		// 	cerr << "Refining alignment, step : " << iter_ctr << endl;

//			// For debugging why 1fyvA00 vs 2rirA01 never stops
//			const protein &protein_a         = prm_proteins[0];
//			const protein &protein_b         = prm_proteins[1];
//			const path temp_align_out_file( "temp_align." + lexical_cast<string>(iter_ctr) + ".txt" );
//			ofstream temp_align_ofstream;
//			open_ofstream( temp_align_ofstream, temp_align_out_file );
//			output_alignment_to_cath_ssap_legacy_format( temp_align_ofstream, next_alignment, protein_a, protein_b );
//			temp_align_ofstream.close();
		// }

		if ( next_matches_prev ) {
			break;
		}

		// ++iter_ctr;
	}
//	const protein &protein_a         = prm_proteins[0];
//	const protein &protein_b         = prm_proteins[1];
//	score_alignment( residue_scorer(), curr_alignment, prm_proteins );
//	output_alignment_to_cath_ssap_legacy_format( cout, curr_alignment, protein_a, protein_b );

	return curr_alignment;
}


/// \brief TODOCUMENT
alignment alignment_refiner::iterate_join(const alignment    &prm_alignment,   ///< TODOCUMENT
                                          const protein_list &prm_proteins,    ///< TODOCUMENT
                                          const gap_penalty  &prm_gap_penalty, ///< TODOCUMENT
                                          const size_vec     &prm_group        ///< TODOCUMENT
                                          ) {
	return iterate_join(
		prm_alignment,
		prm_proteins,
		view_cache_list( prm_proteins ),
		prm_gap_penalty,
		prm_group
	);
}

/// \brief TODOCUMENT
alignment alignment_refiner::iterate_join(const alignment       &prm_alignment,       ///< TODOCUMENT
                                          const protein_list    &prm_proteins,        ///< TODOCUMENT
                                          const view_cache_list &prm_view_cache_list, ///< TODOCUMENT
                                          const gap_penalty     &prm_gap_penalty,     ///< TODOCUMENT
                                          const size_vec        &prm_group            ///< TODOCUMENT
                                          ) {
	return iterate_step_for_alignment_split_list(
		prm_alignment,
		prm_proteins,
		prm_view_cache_list,
		prm_gap_penalty,
		make_list_of_alignment_split( prm_alignment, prm_group )
	).second;
}
