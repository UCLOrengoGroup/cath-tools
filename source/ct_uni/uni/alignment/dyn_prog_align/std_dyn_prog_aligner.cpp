/// \file
/// \brief The std_dyn_prog_aligner class definitions

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

#include "std_dyn_prog_aligner.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include "alignment/alignment.hpp"
#include "alignment/dyn_prog_align/detail/matrix_plotter/gnuplot_matrix_plotter.hpp" // ***** TEMPORARY *****
#include "alignment/dyn_prog_align/detail/matrix_plotter/matrix_plot.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/difference.hpp"
#include "common/type_aliases.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::align::gap;
using namespace cath::common;
using namespace std;

using boost::adaptors::reversed;
using boost::numeric_cast;

/// \brief A standard do_clone method.
unique_ptr<dyn_prog_aligner> std_dyn_prog_aligner::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Align two sequences of items to maximise their match scores
///
/// Overview
/// --------
///
/// This is an implementation of the standard NW algorithm (although it uses the O(n^2) refinement
/// rather than the O(n^3) original).
///
/// This uses a dynamic-programming approach - it iteratively solves the optimum alignment from particular points to the end
/// by reusing previously calculated results that start from later points.
///
/// This handles gap-penalties (a cost that is subtracted from an alignment's score each time it opens a gap)
/// but note that, in-line with the previous SSAP code, it repeatedly applies the gap penalty for every step in a gap.
///
/// This also allows the algorithm to be windowed so that it doesn't explore the areas that result in an alignment
/// with a large offset.
///
/// The setup
/// ---------
///
/// Remember that this is a two-dimensional "fence-post problem" (look it up to get the basic idea)
/// where the scores are fence panels and points through which the path traces are fence posts.
///
/// To illustrate with the following diagram:
///  * the numbers in the squares represent entries in the score matrix
///    (which are numbered from (0, 0) to (n - 1, m - 1) )
///  * the corners of the squares represent the points through which the alignment
///    can pass (horizontally, vertically or down-one-and-right-one diagonally)
///  * the hashes represent entries in the return_path matrix
///    (which are numbered from (0, 0) to (n - 1, m - 1) )
///
///       .       0     1     2    ---------->   m-3   m-2   m-1
///
///            #-----#-----#-----#             #-----#-----#-----+
///       0    |  5  |  1  |  1  | . . . . . . |  1  |  0  |  2  |
///            #-----#-----#-----#             #-----#-----#-----+
///       1    |  3  |  6  |  2  | . . . . . . |  0  |  0  |  2  |
///            #-----#-----#-----#             #-----#-----#-----+
///       2    |  8  |  3  |  4  | . . . . . . |  0  |  1  |  0  |
///            #-----#-----#-----#             #-----#-----#-----+
///       |       .     .     .    .              .     .     .
///       |       .     .     .      .            .     .     .
///       |       .     .     .        .          .     .     .
///       |       .     .     .          .        .     .     .
///       |       .     .     .            .      .     .     .
///       V       .     .     .              .    .     .     .
///            #-----#-----#-----#             #-----#-----#-----+
///      n-3   |  1  |  0  |  2  | . . . . . . |  6  |  9  |  5  |
///            #-----#-----#-----#             #-----#-----#-----+
///      n-2   |  2  |  1  |  0  | . . . . . . |  3  |  8  |  7  |
///            #-----#-----#-----#             #-----#-----#-----+
///      n-1   |  0  |  1  |  0  | . . . . . . |  6  |  8  |  7  |
///            +-----+-----+-----+             +-----+-----+-----+
///
/// A valid path runs from the very top left to the very bottom right in a series of straight
/// steps between squares' corners. The permitted steps are: down-one, right-one or diagonally
/// down-one-and-right-one. In practice, this implementation should never produce a path with right angles.
///
/// The return_path matrix describes the optimum path that has been found from any
/// given position to the bottom right. The algorithm consists of two steps:
///  1. populating the return_path matrix from bottom right to top left
///  1. tracing back along the optimum path to the top-left
///
/// Populating the return_path matrix
/// ---------------------------------
///
/// The code sweeps from bottom-right to top-left, working up each column before moving
/// one column to the left. Since the matrix may be windowed, the sweep up each column
/// may only cover a sub-strip of the cells.
///
/// For each point in the return_path matrix (ie each hash in the above figure), the optimum
/// path to the end can be calculated by considering three paths:
///  * move one step down-and-right and then take the best (pre-calculated) path from there (incurring no gap penalty and potentially adding extra score).
///  * move one step down           and then take the best path from there (incurring one gap penalty and adding no extra score)
///  * move one step right          and then take the best path from there (incurring one gap penalty and adding no extra score)
///
/// If the scores are equally good, then the code will prefer the first (diagonal) move and then prefer whichever it has been instructed to prefer.
///
/// ABEDE vs ACDED with gap penalty 1:
///    A  B  E  D  E
/// A  1  0  0  0  0
/// C  0  0  0  0  0
/// D  0  0  0  1  0
/// E  0  0  1  0  1
/// D  0  0  0  1  0
///
/// ABEDE  has score 2
///  ACDED
///
///  ABEDE has score 2
/// ACDED
///
/// A BEDE has score 2
/// ACDED
///
/// ABEDE  has score 2
/// A CDED
///
/// Sometimes there is a choice between multiple different return paths that achieve the same score,
/// In these cases the code uses the following criteria to attempt to distinguish (in descending order):
///  * prefer a diagonal move
///  * prefer a move toward the diagonal passing through the top-left corner
///  * prefer a move toward the diagonal passing through the bottom-right corner
///  * fall back on a parameter that instructs which way to go

score_alignment_pair std_dyn_prog_aligner::do_align(const dyn_prog_score_source &prm_scorer,          ///< TODOCUMENT
                                                    const gap_penalty           &prm_gap_penalty,     ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
                                                    const size_type             &/*prm_window_width*/ ///< TODOCUMENT
                                                    ) const {
	const size_t length_a     = prm_scorer.get_length_a();
	const size_t length_b     = prm_scorer.get_length_b();
	const size_t window_width = get_window_width_for_full_matrix(length_a, length_b);

	the_return_path.reset(        length_a, length_b, window_width );
	the_accumulated_scores.reset( length_a, length_b, window_width );

	for (const size_t &index_a : indices( length_a ) | reversed ) {
		for (const size_t &index_b : indices( length_b ) | reversed ) {

			const path_step_score_map score_of_path_step = get_total_scores_of_path_steps_from_point(
				the_accumulated_scores,
				the_return_path,
				prm_gap_penalty,
				prm_scorer,
				index_a,
				index_b
			);

			const path_step the_chosen_path = choose_path_step(
				score_of_path_step,
				index_a,
				index_b,
				length_a,
				length_b
			);
			const score_type the_chosen_score = score_of_path_step.at( the_chosen_path );
			the_return_path.set_path_step_towards_end_at_point(    index_a, index_b, the_chosen_path  );
			the_accumulated_scores.set_score_towards_end_at_point( index_a, index_b, the_chosen_score );

//			cerr << "At ";
//			cerr << index_a;
//			cerr << ", ";
//			cerr << index_b;
//			cerr << ", the chosen path_step is ";

//			matrix_plot<gnuplot_matrix_plotter>( prm_scorer, the_return_path, the_accumulated_scores );
		}
	}
//	matrix_plot<gnuplot_matrix_plotter>( "std_dyn_prog_aligner", prm_scorer, the_return_path, the_accumulated_scores );
	const alignment  final_alignment( make_alignment(the_return_path)                             );
	const score_type final_score(     the_accumulated_scores.get_score_towards_end_at_point(0, 0) );
	return make_pair( final_score, final_alignment );
}

/// \brief TODOCUMENT
path_step std_dyn_prog_aligner::choose_path_step(const path_step_score_map &prm_score_of_path, ///< TODOCUMENT
                                                 const size_t              &prm_index_a,       ///< TODOCUMENT
                                                 const size_t              &prm_index_b,       ///< TODOCUMENT
                                                 const size_t              &prm_length_a,      ///< TODOCUMENT
                                                 const size_t              &prm_length_b       ///< TODOCUMENT
                                                 ) const {
	// Grab the maximum score that's achieved by any of the possible path_steps
	const score_type max_score = max_path_step_score(prm_score_of_path);

//	cerr << "At ";
//	cerr << prm_index_a;
//	cerr << ", ";
//	cerr << prm_index_b;
//	cerr << " [";
//	for (const path_step_score_pair &path_and_score : prm_score_of_path) {
//		cerr << " ";
//		cerr << path_and_score.first;
//		cerr << ":";
//		cerr << path_and_score.second;
//	}
//	cerr << " ] max_score:";
//	cerr << max_score;
//	cerr << endl;

//	const score_type align_pair_score         = prm_score_of_path.at( path_step::ALIGN_PAIR         );
	const score_type insert_into_first_score  = prm_score_of_path.at( path_step::INSERT_INTO_FIRST  );
	const score_type insert_into_second_score = prm_score_of_path.at( path_step::INSERT_INTO_SECOND );

	// If only path_step::ALIGN_PAIR achieves the maximum score then prefer that
	if ( insert_into_first_score != max_score && insert_into_second_score != max_score ) {
//		cerr << "Returning path_step::ALIGN_PAIR" << endl;
		return path_step::ALIGN_PAIR;
	}
	// Else, if path_step::INSERT_INTO_SECOND doesn't achieve the maximum score but path_step::INSERT_INTO_FIRST does, then choose that
	if ( insert_into_first_score == max_score && insert_into_second_score != max_score ) {
//		cerr << "Returning path_step::INSERT_INTO_FIRST" << endl;
		return path_step::INSERT_INTO_FIRST;
	}
	// Else, if path_step::INSERT_INTO_FIRST doesn't achieve the maximum score but path_step::INSERT_INTO_SECOND does, then choose that
	if ( insert_into_first_score != max_score && insert_into_second_score == max_score ) {
//		cerr << "Returning path_step::INSERT_INTO_SECOND" << endl;
		return path_step::INSERT_INTO_SECOND;
	}

	// Otherwise, path_step::INSERT_INTO_FIRST and path_step::INSERT_INTO_SECOND achieve the same maximum score,
	//
	// The remaining code must choose between these two options, preferably in a way
	// that's symmetric in the sense that the same (mirrored) alignment will
	// be generated by flipping the order of the two inputs.

	const map<path_step, size_size_pair> after_indices_by_step = indices_of_point_by_path_step(
		prm_index_a,
		prm_index_b
	);

	// Prefer the step that moves closer to the line between the start and end
	const double first_to_second_ratio =   numeric_cast<double>(prm_length_a)
	                                     / numeric_cast<double>(prm_length_b);
	const size_size_pair indices_after_insert_into_first  = after_indices_by_step.at( path_step::INSERT_INTO_FIRST  );
	const size_size_pair indices_after_insert_into_second = after_indices_by_step.at( path_step::INSERT_INTO_SECOND );
	const double ratio_after_insert_into_first  =   numeric_cast<double>( indices_after_insert_into_first.first   )
	                                              / numeric_cast<double>( indices_after_insert_into_first.second  );
	const double ratio_after_insert_into_second =   numeric_cast<double>( indices_after_insert_into_second.first  )
		                                          / numeric_cast<double>( indices_after_insert_into_second.second );
	const double radio_diff_after_insert_into_first  = difference( ratio_after_insert_into_first,  first_to_second_ratio );
	const double radio_diff_after_insert_into_second = difference( ratio_after_insert_into_second, first_to_second_ratio );

	if ( radio_diff_after_insert_into_first  < radio_diff_after_insert_into_second ) {
		return path_step::INSERT_INTO_FIRST;
	}
	if ( radio_diff_after_insert_into_second < radio_diff_after_insert_into_first  ) {
		return path_step::INSERT_INTO_SECOND;
	}

	const double radio_offset_after_insert_into_first  = difference( ratio_after_insert_into_first,  1.0 );
	const double radio_offset_after_insert_into_second = difference( ratio_after_insert_into_second, 1.0 );
	if ( radio_offset_after_insert_into_first < radio_offset_after_insert_into_second ) {
		return path_step::INSERT_INTO_FIRST;
	}
	if ( radio_diff_after_insert_into_second  < radio_offset_after_insert_into_first  ) {
		return path_step::INSERT_INTO_SECOND;
	}

	return path_step::INSERT_INTO_FIRST;
}
