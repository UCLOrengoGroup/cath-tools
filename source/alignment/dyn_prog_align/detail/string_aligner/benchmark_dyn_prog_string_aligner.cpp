/// \file
/// \brief The benchmark_dyn_prog_string_aligner class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "benchmark_dyn_prog_string_aligner.h"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm/reverse.hpp>

#include "alignment/gap/gap_penalty.h"
#include "exception/invalid_argument_exception.h"
#include "exception/not_implemented_exception.h"

#include <algorithm>
#include <string>

using namespace cath;
using namespace cath::align::gap;
using namespace cath::align::detail;
using namespace cath::common;
using namespace std;

using boost::adaptors::reverse;
using boost::numeric_cast;

/// \brief TODOCUMENT
void benchmark_dyn_prog_string_aligner::make_ends_spaces(string &arg_string ///< TODOCUMENT
                                                         ) {
	for (char &letter : arg_string) {
		if (letter != '-') {
			break;
		}
		letter = ' ';
	}
	for (char &letter : reverse( arg_string) ) {
		if (letter != '-') {
			break;
		}
		letter = ' ';
	}
}

/// \brief TODOCUMENT
///
/// Carry out the alignment of two strings with a specified gap penalty
/// The wording of the comments and variable names etc assumes that sequence1 is written from left to right
/// along the top and sequence2 is written down the left-hand side
///
/// The algorithm scans down each column and then moves to the next column (to the right).
///
/// To match this, the elements of the matrix count up in the same order.
/// Unfortunately, this means that the matrix doesn't use standard matrix notation.
///
/// matrix             - Score matrix (including gap penalties)
/// noOfMatchesMatrix  - Matrix of the number of matches that have been achieved (no gap penalties)
///
/// To record the path taken to make tracing back easier, the following two matrices are used:
/// prevSeq1PosnMatrix - Stores the seq1 position of the previous entry in the path to this entry
/// prevSeq2PosnMatrix - Stores the seq2 position of the previous entry in the path to this entry
///
/// To improve the speed of scanning for the best previous score for a given element, the 4 following values are used:
/// bestColVal         - The best value that has been found so far in the previous column
/// bestColIndex       - The index of the entry containing the best value found so far in the previous column
/// bestRowVal         - An array containing the best values found so far in each row
/// bestRowIndex       - An array containing the indices of the best values found so far in each row
///
/// \todo With gap penalty of 0, this currently aligns DECECC and CAD
///                : CECC
///                : C-AD
///
str_str_pair benchmark_dyn_prog_string_aligner::do_align(const string      &arg_string_a,   ///< TODOCUMENT
                                                         const string      &arg_string_b,   ///< TODOCUMENT
                                                         const gap_penalty &arg_gap_penalty ///< The gap penalty to be applied
                                                         ) const {
	// Check that the open and extend gap penalties are equal and then grab one of them
	if ( arg_gap_penalty.get_open_gap_penalty() != arg_gap_penalty.get_extend_gap_penalty() ) {
		BOOST_THROW_EXCEPTION(not_implemented_exception("benchmark_dyn_prog_string_aligner unable to handle non-equal open and extend gap_penalty values"));
	}
	const score_type   common_gap_penalty = arg_gap_penalty.get_open_gap_penalty();

	const size_t       length_a   = arg_string_a.length();
	const size_t       length_b   = arg_string_b.length();
	const size_t       max_length = max(length_a, length_b);
	vector<score_type> accumulation_matrix(    max_length * max_length );
	vector<ptrdiff_t>  prev_seq_1_posn_matrix( max_length * max_length );
	vector<ptrdiff_t>  prev_seq_2_posn_matrix( max_length * max_length );
	vector<score_type> best_row_value( max_length );
	vector<ptrdiff_t>  best_row_index( max_length );
	score_type         top_score      = 0;
	size_t             best_seq_1_end = 0;
	size_t             best_seq_2_end = 0;

	// Check that the sequences aren't too long
	if (length_a > max_length || length_b > max_length) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"ERROR: align() has been called with a sequence that is longer than the dim of the dynProg"
		));
		exit(1);
	}

	// Loop over all the residues in the first sequence
	for (size_t seq_1_ctr = 0; seq_1_ctr < length_a; ++seq_1_ctr) {
		// If we are on the first residue (index 0), then process separately
		if (seq_1_ctr == 0) {
			for (size_t seq_2_ctr = 0; seq_2_ctr < length_b; ++seq_2_ctr) {
				best_row_index[seq_2_ctr] = 0;

				// If residues seq1Ctr and seq2Ctr are a good match then set matrix entries to 1
				if (arg_string_a[seq_1_ctr] == arg_string_b[seq_2_ctr] && arg_string_a[seq_1_ctr] != 'X') {
					accumulation_matrix[seq_2_ctr] = 1;
					best_row_value     [seq_2_ctr] = 1;
					top_score = 1;
					best_seq_1_end = seq_1_ctr;
					best_seq_2_end = seq_2_ctr;
					prev_seq_1_posn_matrix[seq_2_ctr] = 0;
					prev_seq_2_posn_matrix[seq_2_ctr] = 0;
				}
				else {
					// Else if residues seq1Ctr and seq2Ctr aren't a good match then set matrix entries to 0
					accumulation_matrix   [seq_2_ctr] = 0;
					best_row_value        [seq_2_ctr] = 0;
					prev_seq_1_posn_matrix[seq_2_ctr] = 0;
					prev_seq_2_posn_matrix[seq_2_ctr] = 0;
				}
			}
		}

		// Else if we aren't on the first residue (index > 0), then process normally
		else {
			score_type best_col_value = 0;
			ptrdiff_t  best_col_index = 0;
			score_type prev_value     = 0;

			// Loop over the residues in the second sequence
			for (size_t seq_2_ctr = 0; seq_2_ctr < length_b; ++seq_2_ctr) {
				ptrdiff_t prev_seq_1_posn = numeric_cast<ptrdiff_t>( seq_1_ctr ) - 1;
				if (prev_seq_1_posn < 0) {
					prev_seq_1_posn = 0;
				}
				ptrdiff_t prev_seq_2_posn = numeric_cast<ptrdiff_t>( seq_2_ctr ) - 1;
				if (prev_seq_2_posn < 0) {
					prev_seq_2_posn = 0;
				}

				// If not at the the top of a column and if the element to the top left of this element
				// is better than thisVal (may not be due to gap penalties) then start with that as our best find
				score_type this_value = 0;
				if (seq_2_ctr > 0) {
					const score_type score_up_left = accumulation_matrix[ max_length * (seq_1_ctr - 1) + seq_2_ctr - 1 ];
					this_value = max(this_value, score_up_left);
				}

				// If not at the top of a column and something better has been discovered on the preceding row
				// then the current values to point to that entry
				if (seq_2_ctr > 0) {
					const score_type best_val_on_prec_row = best_row_value[ seq_2_ctr - 1 ];
					if (best_val_on_prec_row - common_gap_penalty > this_value) {
						//if (bestRowVal[seq2Ctr-1] - gap_pen > thisVal || noOfMatchesMatrix[dim * (bestRowIndex[seq2Ctr-1]) + seq2Ctr-1] < thisMatchCount) {
						this_value      = best_val_on_prec_row - common_gap_penalty;
						prev_seq_1_posn = best_row_index[seq_2_ctr - 1];
						prev_seq_2_posn = numeric_cast<ptrdiff_t>( seq_2_ctr ) - 1;
						//}
					}
				}

				// If something better has been discovered on the preceding column
				// then the current values to point to that entry
				if (best_col_value - common_gap_penalty > this_value) {
					//if (bestColVal - gap_pen > thisVal || noOfMatchesMatrix[dim * (seq1Ctr-1) + bestColIndex] < thisMatchCount) {
						this_value      = best_col_value - common_gap_penalty;
						prev_seq_1_posn = numeric_cast<ptrdiff_t>( seq_1_ctr ) - 1;
						prev_seq_2_posn = best_col_index;
					//}
				}

				// This cell represents aligning argString1[seq1Ctr] and argString2[seq2Ctr] so
				// if they are a good match, award another point
				if (arg_string_a[seq_1_ctr] == arg_string_b[seq_2_ctr] && arg_string_a[seq_1_ctr] != 'X') {
					++this_value;
				}

				// If this is the best score we've seen, then we should stop here unless we see something better later
				// so record the position and score
				if (this_value >= top_score) {
					top_score      = this_value;
					best_seq_1_end = seq_1_ctr;
					best_seq_2_end = seq_2_ctr;
				}

				// Before moving on to the next cell, update the bestRowVal and the relevant entry in bestCol

				// Update the best column index and value
				score_type score_left = accumulation_matrix[max_length * (seq_1_ctr-1) + seq_2_ctr];
				if (score_left > best_col_value) {
					best_col_value = score_left;
					best_col_index = numeric_cast<ptrdiff_t>( seq_2_ctr );
				}

				// Update the best row index and value for the previous row
				if (seq_2_ctr) {
					if (prev_value > best_row_value[seq_2_ctr - 1]) {
						best_row_value[ seq_2_ctr - 1 ] = prev_value;
						best_row_index[ seq_2_ctr - 1 ] = numeric_cast<ptrdiff_t>( seq_1_ctr );
					}
				}

				// Now we've finished for this element so record the value of the score in the score matrix
				accumulation_matrix   [ max_length * seq_1_ctr + seq_2_ctr ] = this_value;
				prev_seq_1_posn_matrix[ max_length * seq_1_ctr + seq_2_ctr ] = prev_seq_1_posn;
				prev_seq_2_posn_matrix[ max_length * seq_1_ctr + seq_2_ctr ] = prev_seq_2_posn;

				// Store the value as the previous value so that after the next iteration, we can update
				// the best row index and value for this row

				prev_value = this_value;
			}
		}
	}

	// Traceback through the matrix to find the alignment
	string aligned_1;
	string aligned_2;

	// Construct the tails first
	ptrdiff_t seq_1_tail_length = numeric_cast<ptrdiff_t>( length_a ) - numeric_cast<ptrdiff_t>( best_seq_1_end );
	ptrdiff_t seq_2_tail_length = numeric_cast<ptrdiff_t>( length_b ) - numeric_cast<ptrdiff_t>( best_seq_2_end );

	if (seq_1_tail_length > seq_2_tail_length) {
		aligned_2 += string( numeric_cast<size_t>( seq_1_tail_length - seq_2_tail_length ), '-');
	}
	else if (seq_2_tail_length > seq_1_tail_length) {
		aligned_1 += string( numeric_cast<size_t>( seq_2_tail_length - seq_1_tail_length ), '-');
	}
	for (ptrdiff_t tail_1_ctr = numeric_cast<ptrdiff_t>(length_a) - 1; tail_1_ctr >= numeric_cast<ptrdiff_t>( best_seq_1_end ); --tail_1_ctr) {
		aligned_1 += arg_string_a[ numeric_cast<size_t>( tail_1_ctr ) ];
	}
	for (ptrdiff_t tail_2_ctr = numeric_cast<ptrdiff_t>(length_b) - 1; tail_2_ctr >= numeric_cast<ptrdiff_t>( best_seq_2_end ); --tail_2_ctr) {
		aligned_2 += arg_string_b[ numeric_cast<size_t>( tail_2_ctr ) ];
	}

	// Now work back through the matrix, appending to aligned1 and aligned2
	ptrdiff_t seq_1_ctr = numeric_cast<ptrdiff_t>( best_seq_1_end );
	ptrdiff_t seq_2_ctr = numeric_cast<ptrdiff_t>( best_seq_2_end );
	while (seq_1_ctr > 0 && seq_2_ctr > 0) {
		ptrdiff_t new_seq_1_ctr = prev_seq_1_posn_matrix[ max_length * numeric_cast<size_t>( seq_1_ctr ) + numeric_cast<size_t>( seq_2_ctr ) ];
		ptrdiff_t new_seq_2_ctr = prev_seq_2_posn_matrix[ max_length * numeric_cast<size_t>( seq_1_ctr ) + numeric_cast<size_t>( seq_2_ctr ) ];
		//cerr << "(" << seq1Ctr << "," << seq2Ctr << ") -> (" << newSeq1Ctr << "," << newSeq2Ctr << ")" << endl;
		if (seq_1_ctr - new_seq_1_ctr > 1) {
			aligned_2 += string( numeric_cast<size_t>( seq_1_ctr - new_seq_1_ctr - 1 ), '-');
			for (ptrdiff_t ins_1_ctr = seq_1_ctr-1; ins_1_ctr > new_seq_1_ctr; ins_1_ctr--) {
				aligned_1 += arg_string_a[ numeric_cast<size_t>( ins_1_ctr ) ];
			}
		}
		else if (seq_2_ctr - new_seq_2_ctr > 1) {
			aligned_1 += string( numeric_cast<size_t>( seq_2_ctr - new_seq_2_ctr - 1 ), '-');
			for (ptrdiff_t ins_2_ctr = seq_2_ctr-1; ins_2_ctr > new_seq_2_ctr; ins_2_ctr--) {
				aligned_2 += arg_string_b[ numeric_cast<size_t>( ins_2_ctr ) ];
			}
		}
		aligned_2 += arg_string_b[ numeric_cast<size_t>( new_seq_2_ctr ) ];
		aligned_1 += arg_string_a[ numeric_cast<size_t>( new_seq_1_ctr ) ];
		seq_1_ctr = new_seq_1_ctr;
		seq_2_ctr = new_seq_2_ctr;
	}

	// Finally construct the heads
	for (ptrdiff_t head1Ctr = seq_1_ctr-1; head1Ctr >= 0; --head1Ctr) {
		aligned_1 += arg_string_a[ numeric_cast<size_t>( head1Ctr ) ];
	}
	for (ptrdiff_t head2Ctr = seq_2_ctr-1; head2Ctr >= 0; --head2Ctr) {
		aligned_2 += arg_string_b[ numeric_cast<size_t>( head2Ctr ) ];
	}
	if (seq_1_ctr > seq_2_ctr) {
		aligned_2 += string( numeric_cast<size_t>( seq_1_ctr - seq_2_ctr ), '-' );
	}
	else if (seq_2_ctr > seq_1_ctr) {
		aligned_1 += string( numeric_cast<size_t>( seq_2_ctr - seq_1_ctr ), '-' );
	}

	// Now we have traced-back, reverse the strings
	reverse( aligned_1 );
	reverse( aligned_2 );

	make_ends_spaces(aligned_1);
	make_ends_spaces(aligned_2);

	return make_pair(aligned_1, aligned_2);
}

