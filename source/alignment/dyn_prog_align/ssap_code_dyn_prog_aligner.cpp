/// \file
/// \brief The ssap_code_dyn_prog_aligner class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 1989, Orengo Group, University College London
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

#include "ssap_code_dyn_prog_aligner.h"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.h"
#include "alignment/gap/gap_penalty.h"
#include "alignment/pair_alignment.h"
#include "common/clone/make_uptr_clone.h"
#include "common/debug_numeric_cast.h"
#include "common/type_aliases.h"
#include "exception/not_implemented_exception.h"
#include "exception/out_of_range_exception.h"
#include "ssap/windowed_matrix.h"

#include <cassert> // **** TEMPORARY *****

using namespace cath;
using namespace cath::align;
using namespace cath::align::gap;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;
using boost::numeric_cast;

/// \brief A standard do_clone method.
unique_ptr<dyn_prog_aligner> ssap_code_dyn_prog_aligner::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Determines the best pathways through the upper and lower score matrices
///        using the Needleman and Wunsch dynamic programming algorithm
///
/// This finds the best pathway through the matrix and returns that alignment along with the score.
///
/// \todo Continue sorting out this nightmare of a subroutine, which was a huge violator of the
///       Dependency Inversion Principle (http://en.wikipedia.org/wiki/Dependency_inversion_principle)
///       but which has gotten a bit simpler since populate_upper_score_matrix() was separated out.
///
/// Notes on trying to figure this out (~September 2013)
/// ====================================================
///
/// Overview
/// --------
///
/// This currently operates on the upper and lower matrices.
/// It doesn't store its own score matrix, it just retrieves values as it needs them.
///
/// Organisation of dynamic programming (DP)
/// ----------------------------------------
///
/// The dynamic programming (DP) code nomenclature considers the matrix of a versus b scores
/// to be laid out like this:
///
///             ---b-->
///          + + + + + + +
///       |  + + + + + + +
///       a  + + + + + + +
///       |  + + + + + + +
///       V  + + + + + + +
///          + + + + + + +
///
/// Both dimensions (currently) use an offset of 1, ie they're labelled from 1 to length (inclusive).
///
/// Note that the actual matrix is indexed by b and then a (which I think someone has done
/// in the past to make existing indexing errors less catastrophic) but just ignore this when
/// considering rows/columns.
/// \todo Rectify this as part of moving to use a windowed matrix class
///
/// The DP code sweeps from bottom-right to top-left, working up each column before moving
/// one column to the left. Since the matrix is typically windowed, the sweep up each column
/// may only cover a sub-strip of the cells.
///
/// Still fairly opaque:
/// --------------------
///
/// It is not clear how the following variables are used.
///  - enter
///  - rat
///
///
///
ssap_code_dyn_prog_aligner::size_size_int_int_score_tuple ssap_code_dyn_prog_aligner::score_matrix(const dyn_prog_score_source &arg_scorer,       ///< TODOCUMENT
                                                                                                   const score_type            &arg_gap_penalty,  ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
                                                                                                   const size_type             &arg_window_width, ///< TODOCUMENT
                                                                                                   int_vec_vec                 &arg_path_matrix   ///< TODOCUMENT
                                                                                                   ) {
	const size_t     &length_a           = arg_scorer.get_length_a();
	const size_t     &length_b           = arg_scorer.get_length_b();

	const score_type  VERY_POOR_SCORE    = numeric_limits<score_type>::min() / 10;
	score_type        best_score         = VERY_POOR_SCORE;
	size_t            flip_flop_current  = 0;
	size_t            flip_flop_previous = 1;
	size_t            mat_a              = 0;
	size_t            mat_b              = 0;
	int               final_path_a       = 0;
	int               final_path_b       = 0;
	int               enter              = 0;
	bool              edge_set           = false;

	// The best scores...???
	/// \todo Are the +2s necessary?
	static score_vec best_scores_in_column;
	best_scores_in_column.assign( arg_window_width + 2, 0 );

	// The indices corresponding to the best scores...???
	/// \todo Are the +2s necessary?
	static size_vec indices_of_best_scores_in_column;
	indices_of_best_scores_in_column.assign( arg_window_width + 2, 0 );

	// Matrix to store row scores in a flip-flop fashion (ie two sets of values: one active; one inactive)
	/// \todo Are the +2s necessary?
	static score_vec_vec row_scores_flipflop_matrix;
	row_scores_flipflop_matrix.assign( 2, score_vec( arg_window_width + 2, VERY_POOR_SCORE ) );

	// Initialise various variable for the right-most column
	for (size_t a_dest_to_index = 0; a_dest_to_index <= arg_window_width + 1; ++a_dest_to_index) {
		row_scores_flipflop_matrix[ 0 ][ a_dest_to_index ] = VERY_POOR_SCORE;
		row_scores_flipflop_matrix[ 1 ][ a_dest_to_index ] = VERY_POOR_SCORE;
		arg_path_matrix[ a_dest_to_index ][ length_b ]     = 0;
		best_scores_in_column[ a_dest_to_index ]           = VERY_POOR_SCORE;
	}

	// Compare each element in protein B with each element in protein A within window
	for (size_t ctr_b = length_b; ctr_b > 0; --ctr_b) {
		score_type best_row_score = VERY_POOR_SCORE;
		size_t     best_row_index = 0;

		// Set pointer to element in protein B
		const size_t window_start = get_window_start_a_for_b__offset_1( length_a, length_b, arg_window_width, ctr_b );
		const size_t window_stop  = get_window_stop_a_for_b__offset_1 ( length_a, length_b, arg_window_width, ctr_b );

		// Set window range, if comparing residues, window is set by selected residues
		if ( --enter < 0 ) {
			enter = debug_numeric_cast<int>(arg_window_width) - 1;
		}

		if ( window_stop == length_a ) {
			arg_path_matrix[length_a - window_start][ctr_b] = 0;
		}

		for (size_t ctr_a = window_stop; ctr_a >= window_start; --ctr_a ) {
			const int a_matrix_idx = get_window_matrix_a_index__offset_1(length_a, length_b, arg_window_width, ctr_a, ctr_b);
			int       rat          = enter + a_matrix_idx;

//			cerr << "Getting score from " << ctr_a << " (os1) and " << ctr_b << " (os1) : " << get_score__offset_1(arg_scorer, ctr_a, ctr_b) << endl;
			row_scores_flipflop_matrix[flip_flop_current][ numeric_cast<size_t>( a_matrix_idx ) ] = get_score__offset_1( arg_scorer, ctr_a, ctr_b );

			if ( ctr_a == length_a || ctr_b == length_b ) {
				continue;
			}

			// ACCUMULATING MATRIX

			// Adjust rat index and window edge
			if ( rat < 0 ) {
				rat += debug_numeric_cast<int>(arg_window_width);
			}
			if ( rat >= debug_numeric_cast<int>(arg_window_width) ) {
				rat -= debug_numeric_cast<int>(arg_window_width);
			}
			if (ctr_a == window_start && !edge_set) {
				best_scores_in_column           [ numeric_cast<size_t>( rat ) ] = VERY_POOR_SCORE;
				indices_of_best_scores_in_column[ numeric_cast<size_t>( rat ) ] = ctr_b;
				if ( ctr_a == 1 ) {
					edge_set = true;
				}
			}

			// Set diagonal score and maximum row and column scores
			const score_type best_col_score = best_scores_in_column           [ numeric_cast<size_t>( rat ) ];
			const size_t     best_col_index = indices_of_best_scores_in_column[ numeric_cast<size_t>( rat ) ];

			const score_type diag_score     = row_scores_flipflop_matrix[flip_flop_previous][ numeric_cast<size_t>( a_matrix_idx ) ];
			const score_type col_score      = best_col_score - arg_gap_penalty;
			const score_type row_score      = best_row_score - arg_gap_penalty;

			// If diagonal cell score greater than max score from row or column - penalty,
			// accumulate diagonal score
			if ( diag_score >= col_score && diag_score >= row_score ) {
				row_scores_flipflop_matrix[flip_flop_current][ numeric_cast<size_t>( a_matrix_idx ) ] += diag_score;
				arg_path_matrix[ numeric_cast<size_t>( a_matrix_idx ) ][ctr_b]                         = 1;
			}
			// Else if row score is better than column score, accumulate maximum score from row
			else if ( row_score > col_score ) {
				row_scores_flipflop_matrix[flip_flop_current][ numeric_cast<size_t>( a_matrix_idx ) ] +=   row_score;
				arg_path_matrix[ numeric_cast<size_t>( a_matrix_idx ) ][ctr_b]                         =   numeric_cast<int>( best_row_index - ctr_a + 1 );
			}
			// Else accumulate maximum score from column
			else {
				row_scores_flipflop_matrix[flip_flop_current][ numeric_cast<size_t>( a_matrix_idx ) ] +=   col_score;
				arg_path_matrix[ numeric_cast<size_t>( a_matrix_idx ) ][ctr_b]                         = - numeric_cast<int>( best_col_index - ctr_b + 1 );
			}

			// If diagonal score greater than previous maximum for row or column, save
			if (diag_score > best_row_score) {
				best_row_score = diag_score;
				best_row_index = ctr_a;
			}
			if (diag_score > best_col_score) {
				best_scores_in_column[ numeric_cast<size_t>( rat ) ]            = diag_score;
				indices_of_best_scores_in_column[ numeric_cast<size_t>( rat ) ] = ctr_b;
			}

			// Save highest score in matrix and cell coordinates
			if ( row_scores_flipflop_matrix[flip_flop_current][ numeric_cast<size_t>( a_matrix_idx ) ] >= best_score ) {
				best_score   = row_scores_flipflop_matrix[flip_flop_current][ numeric_cast<size_t>( a_matrix_idx ) ];
				final_path_a = a_matrix_idx;
				final_path_b = numeric_cast<int>( ctr_b );
				mat_a        = 1;
				mat_b        = 1;
				if ( ctr_a == 1 ) {
					mat_b = ctr_b;
				}
				if ( ctr_b == 1 ) {
					mat_a = ctr_a;
				}
			}
		}
		swap(flip_flop_current, flip_flop_previous);
	}

	return make_tuple(
		mat_a,
		mat_b,
		final_path_a,
		final_path_b,
		best_score
	);
}

/// \brief Build an alignment by tracing back through a score matrix
///
/// This is a wrapper providing access to trace_recursive().
alignment ssap_code_dyn_prog_aligner::traceback(const size_t      &arg_length_a,     ///< The number of entries (residues or secondary structures) in the first  structure
                                                const size_t      &arg_length_b,     ///< The number of entries (residues or secondary structures) in the second structure
                                                const int         &arg_mat_a,        ///< TODOCUMENT
                                                const int         &arg_mat_b,        ///< TODOCUMENT
                                                const int         &arg_path_index_a, ///< The current a index in the path matrix
                                                const int         &arg_path_index_b, ///< The current b index in the path matrix
                                                const int_vec_vec &arg_path_matrix   ///< A matrix indicating the best step back from each point
                                                ) {
	// Construct an empty alignment
	alignment new_alignment(alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT);

//	// The alignment might be up to the sum of the two lengths so reserve that much memory
	new_alignment.reserve(arg_length_a + arg_length_b);

	// Trace back through matrix from highest scoring cell
	traceback_recursive(
		new_alignment,
		arg_mat_a,
		arg_mat_b,
		arg_path_index_a,
		arg_path_index_b,
		0,
		0,
		arg_path_matrix
	);

#ifndef NDEBUG

	/// \todo Investigate why this sometimes overruns, eg for 2h7cB00 vs 1qo0D02
	///
	/// Quite a few other bugs appear when the lengths are very different, as they are
	/// in this case (531 and 46 respectively).
	///
	/// This appears to be an existing problem in the DP/tracing code because adding
	/// a few debug statements into the original SSAP code and then running on 2h7cB00/1qo0D02
	/// produced 60 copies of the line:
	///
	///     After calling trace, posa[aln->length] : 777, lena: 531, posb[aln->length] : 46, lenb: 46
	if (get_last_a_offset_1_position(new_alignment) > arg_length_a) {
//		cerr << endl;
//		for (size_t ctr_a = 0; ctr_a <= length_a; ++ctr_a) {
//			for (size_t ctr_b = 0; ctr_b <= length_b; ++ctr_b) {
//				cerr << path_matrix[ctr_a][ctr_b] << "\t";
//			}
//			cerr << endl;
//		}
//		cerr << endl;
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"After tracing, the last a_position is "
			+ lexical_cast<string>(get_last_a_offset_1_position(new_alignment))
			+ " which is past the end of structure a, which is of length "
			+ lexical_cast<string>(arg_length_a)
		));
	}
	if (get_last_b_offset_1_position( new_alignment ) > arg_length_b) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"After tracing, the last b_position is "
			+ lexical_cast<string>(get_last_b_offset_1_position( new_alignment ) )
			+ " which is past the end of structure b, which is of length "
			+ lexical_cast<string>(arg_length_b)
		));
	}
#endif

	// Add any overhangs
	// (trace_recursive() will have built the alignment to the end of one of the lists,
	// but there may be more on the other)
	if ( has_last_b_position( new_alignment ) && get_last_b_offset_1_position( new_alignment ) == arg_length_b ) {
		while ( ! has_last_a_position( new_alignment ) || get_last_a_offset_1_position( new_alignment ) < arg_length_a) {
			append_position_a_offset_1( new_alignment, get_last_a_offset_1_position( new_alignment ) + 1 );
		}
	}
	else if ( has_last_a_position( new_alignment ) && get_last_a_offset_1_position( new_alignment ) == arg_length_a ) {
		while ( ! has_last_b_position( new_alignment ) || get_last_b_offset_1_position( new_alignment ) < arg_length_b) {
			append_position_b_offset_1( new_alignment, get_last_b_offset_1_position( new_alignment ) + 1 );
		}
	}

	// Return the new alignment
	return new_alignment;
}


/// \brief Recursively trace back through score matrix to find an optimum pathway
///
/// This should be accessed via the more user-friendly trace() subroutine.
///
/// From what I can gather, pat is a matrix containing the path back from each point
/// where each entry is:
///  1 : go back ???
void ssap_code_dyn_prog_aligner::traceback_recursive(alignment         &arg_alignment,      ///< The alignment being recursively constructed.
                                                     const int         &arg_mat_a,          ///< TODOCUMENT
                                                     const int         &arg_mat_b,          ///< TODOCUMENT
                                                     const int         &arg_path_index_a,   ///< The current i index in the path matrix
                                                     const int         &arg_path_index_b,   ///< The current k index in the path matrix
                                                     const int         &arg_previous_mat_a, ///< The mat_i value in the previous trace() call in the recursive stack
                                                     const int         &arg_previous_mat_b, ///< The mat_j value in the previous trace() call in the recursive stack
                                                     const int_vec_vec &arg_path_matrix     ///< A matrix indicating the best step back from each point
                                                     ) {
	//
	const int path_entry = arg_path_matrix[ numeric_cast<size_t>( arg_path_index_a ) ]
	                                      [ numeric_cast<size_t>( arg_path_index_b ) ];

	// Insertions in sequence A from arg_previous_mat_i + 1 to arg_mat_i - 1 (inclusive)
	for (int ctr_a = arg_previous_mat_a + 1; ctr_a < arg_mat_a; ++ctr_a) {
		append_position_a_offset_1( arg_alignment, numeric_cast<size_t>( ctr_a ) );
	}

	// Insertions in sequence B  from arg_previous_mat_j + 1 to arg_mat_j - 1 (inclusive)
	for (int ctr_b = arg_previous_mat_b + 1; ctr_b < arg_mat_b; ++ctr_b) {
		append_position_b_offset_1( arg_alignment, numeric_cast<size_t>( ctr_b ) );
	}

	// Current aligned position
	append_position_both_offset_1( arg_alignment, numeric_cast<size_t>( arg_mat_a ), numeric_cast<size_t>( arg_mat_b ) );

	//
	if ( !path_entry ) {
		return;
	}

	const int next_mat_a        = arg_mat_a         + (   (path_entry == 1) ?  1
	                                                    : (path_entry <  1) ?  1
	                                                                        :  0 + path_entry );
	const int next_mat_b        = arg_mat_b         + (   (path_entry == 1) ?  1
	                                                    : (path_entry <  1) ?  0 - path_entry
	                                                                        :  1              );
	const int next_path_index_a = arg_path_index_a  + (   (path_entry == 1) ?  0
	                                                    : (path_entry <  1) ?  1 + path_entry
	                                                                        : -1 + path_entry );
	const int next_path_index_b = arg_path_index_b  + (   (path_entry == 1) ?  1
	                                                    : (path_entry <  1) ?  0 - path_entry
	                                                                        :  1              );
	traceback_recursive(
		arg_alignment,
		next_mat_a,
		next_mat_b,
		next_path_index_a,
		next_path_index_b,
		arg_mat_a,
		arg_mat_b,
		arg_path_matrix
	);
}


/// \brief Determines the best pathways through the upper and lower score matrices
///        using the Needleman and Wunsch dynamic programming algorithm
///
/// This finds the best pathway through the matrix and returns that alignment along with the score.
///
/// \todo Continue sorting out this nightmare of a subroutine, which was a huge violator of the
///       Dependency Inversion Principle (http://en.wikipedia.org/wiki/Dependency_inversion_principle)
///       but which has gotten a bit simpler since populate_upper_score_matrix() was separated out.
///
/// Notes on trying to figure this out (~September 2013)
/// ====================================================
///
/// Overview
/// --------
///
/// This currently operates on the upper and lower matrices.
/// It doesn't store its own score matrix, it just retrieves values as it needs them.
///
/// Organisation of dynamic programming (DP)
/// ----------------------------------------
///
/// The dynamic programming (DP) code nomenclature considers the matrix of a versus b scores
/// to be laid out like this:
///
///             ---b-->
///          + + + + + + +
///       |  + + + + + + +
///       a  + + + + + + +
///       |  + + + + + + +
///       V  + + + + + + +
///          + + + + + + +
///
/// Both dimensions (currently) use an offset of 1, ie they're labelled from 1 to length (inclusive).
///
/// Note that the actual matrix is indexed by b and then a (which I think someone has done
/// in the past to make existing indexing errors less catastrophic) but just ignore this when
/// considering rows/columns.
/// \todo Rectify this as part of moving to use a windowed matrix class
///
/// The DP code sweeps from bottom-right to top-left, working up each column before moving
/// one column to the left. Since the matrix is typically windowed, the sweep up each column
/// may only cover a sub-strip of the cells.
///
/// Still fairly opaque:
/// --------------------
///
/// It is not clear how the following variables are used.
///  - enter
///  - rat
///
///
///
score_alignment_pair ssap_code_dyn_prog_aligner::do_align(const dyn_prog_score_source &arg_scorer,      ///< TODOCUMENT
                                                          const gap_penalty           &arg_gap_penalty, ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
                                                          const size_type             &arg_window_width ///< TODOCUMENT
                                                          ) const {
	if (arg_gap_penalty.get_extend_gap_penalty() != 0) {
		BOOST_THROW_EXCEPTION(not_implemented_exception("ssap_code_dyn_prog_aligner unable to handle non-zero extend_gap_penalty"));
	}

	const size_t length_a = arg_scorer.get_length_a();
	const size_t length_b = arg_scorer.get_length_b();

	// Matrix to store the first step in the best path from each cell to the bottom right of the matrix
	/// \todo Are the +2s necessary?
	/// \todo Is the +1 necessary?
	static int_vec_vec path_matrix;
	path_matrix.assign( arg_window_width + 2, int_vec( length_b + 1, 0 ) );

	// Score the matrix and hence build up a matrix of the best path back
	const size_size_int_int_score_tuple score_nums = score_matrix(arg_scorer, arg_gap_penalty.get_open_gap_penalty(), arg_window_width, path_matrix);
	const size_t     &mat_a        = get<0>( score_nums );
	const size_t     &mat_b        = get<1>( score_nums );
	const int        &final_path_a = get<2>( score_nums );
	const int        &final_path_b = get<3>( score_nums );
	const score_type &best_score   = get<4>( score_nums );

	// Trace back through matrix from highest scoring cell
	alignment new_alignment = traceback(
		length_a,
		length_b,
		static_cast<int>( mat_a ),
		static_cast<int>( mat_b ),
		final_path_a,
		final_path_b,
		path_matrix
	);

	return make_pair(best_score, new_alignment);
}

