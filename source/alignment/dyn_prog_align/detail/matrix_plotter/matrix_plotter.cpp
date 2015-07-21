/// \file
/// \brief The matrix_plotter class definitions

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

#include "matrix_plotter.h"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/irange.hpp>

#include "alignment/dyn_prog_align/detail/return_path_matrix.h"
#include "alignment/dyn_prog_align/detail/score_accumulation_matrix.h"
#include "alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.h"
#include "common/boost_addenda/range/utility/iterator/cross_itr.h"
#include "common/c++14/cbegin_cend.h"
#include "exception/invalid_argument_exception.h"
#include "ssap/windowed_matrix.h"

#include <limits>

using namespace boost::filesystem;
using namespace cath::align::detail;
using namespace cath::common;
using namespace std;

using boost::irange;
using boost::numeric_cast;

/// \brief TODOCUMENT
void matrix_plotter::check_lengths_match(const size_t &arg_length_a, ///< TODOCUMENT
                                         const size_t &arg_length_b  ///< TODOCUMENT
                                         ) const {
	if (arg_length_a != get_length_a()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot plot entry associated with a mismatching first length"));
	}
	if (arg_length_b != get_length_b()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot plot entry associated with a mismatching second length"));
	}
}

/// \brief TODOCUMENT
void matrix_plotter::check_window_width_matches(const size_t &arg_window_width ///< TODOCUMENT
                                                ) {
	if (arg_window_width != get_window_width()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot plot entry associated with a mismatching window width"));
	}
}

/// \brief Get the length of the first entry (the number of rows in the matrix)
///
/// This is protected rather than private so that concrete derived classes can use it
size_t matrix_plotter::get_length_a() const {
	return length_a;
}

/// \brief Get the length of the second entry (the number of columns in the matrix)
///
/// This is protected rather than private so that concrete derived classes can use it
size_t matrix_plotter::get_length_b() const {
	return length_b;
}

/// \brief Get the width of the window to be used when plotting
///
/// This is protected rather than private so that concrete derived classes can use it
size_t matrix_plotter::get_window_width() const {
	return window_width;
}

/// \brief Ctor for matrix_plotter
matrix_plotter::matrix_plotter(const size_t &arg_length_a, ///< The length of the first  entry (the number of rows    in the matrix)
                               const size_t &arg_length_b  ///< The length of the second entry (the number of columns in the matrix)
                               ) : length_a ( arg_length_a ),
                                   length_b ( arg_length_b ),
                                   simplified_windowless ( true ) {
}

/// \brief Ctor for matrix_plotter
matrix_plotter::matrix_plotter(const size_t &arg_length_a,    ///< The length of the first  entry (the number of rows    in the matrix)
                               const size_t &arg_length_b,    ///< The length of the second entry (the number of columns in the matrix)
                               const size_t &arg_window_width ///< The window width to use when plotting
                               ) : length_a(arg_length_a),
                                   length_b(arg_length_b),
                                   window_width(arg_window_width) {
}

/// \brief Plot the original scores for the matrix
void matrix_plotter::plot_scores(const dyn_prog_score_source &arg_scorer ///< The source of the scores to plot
                                 ) {
	check_lengths_match( arg_scorer.get_length_a(), arg_scorer.get_length_b() );

	const size_t length_a     = get_length_a();
	const size_t length_b     = get_length_b();
	const size_t window_width = get_window_width();

	doub_vec_vec scores( length_a, doub_vec( length_b, 0.0 ) );

	if ( simplified_windowless ) {
		const auto score_of_idx_tpl = [&] (const tuple<size_t, size_t> &x) { return numeric_cast<double>( arg_scorer.get_score( get<0>( x ), get<1>( x ) ) ); };

		const auto index_range_a     = irange( 0_z, length_a );
		const auto index_range_b     = irange( 0_z, length_b );
		const auto index_range_cross = cross( index_range_a, index_range_b );
		const auto the_minmax = minmax_element(
			common::cbegin( index_range_cross ),
			common::cend  ( index_range_cross ),
			[&] (const tuple<size_t, size_t> &x, const tuple<size_t, size_t> &y) {
				return score_of_idx_tpl( x ) < score_of_idx_tpl( y );
			}
		);

		if ( the_minmax.first == common::cend  ( index_range_cross ) || the_minmax.second == common::cend  ( index_range_cross ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find min/max score in elements"));
		}
		const auto min_score         = score_of_idx_tpl( *the_minmax.first  );
		const auto max_score         = score_of_idx_tpl( *the_minmax.second );

		// Copy scores over
		for (const auto &x : index_range_cross) {
			scores[ get<0>( x ) ][ get<1>( x ) ] = score_of_idx_tpl( x );
		}
		do_plot_scores( scores, min_score, max_score );
	}
	else {

		double min_score = numeric_limits<double>::max();
		double max_score = numeric_limits<double>::min();
		for (size_t ctr_b = 0; ctr_b < length_b; ++ctr_b) {

			/// \todo Tidy up this window stuff:
			///        * remove the need to worry about offset 1 here
			///        * provide a function that returns the start/stop as a pair
			///        * move towards abstracting iterating over windowed matrices into an iterator of windowed_matrix<>
			const size_t window_start_a = get_window_start_a_for_b__offset_1(
				length_a,
				length_b,
				window_width,
				ctr_b + 1
			) - 1;
			const size_t window_stop_a = get_window_stop_a_for_b__offset_1(
				length_a,
				length_b,
				window_width,
				ctr_b + 1
			) - 1;
			for (size_t ctr_a = window_start_a; ctr_a <= window_stop_a; ++ctr_a) {
				const double score = numeric_cast<double>( arg_scorer.get_score( ctr_a, ctr_b ) );
				scores[ctr_a][ctr_b] = score;

				// (pairs flipped due to switching from matrix order (row index then column index) to graph order (x then y))
				if ( score > 0 ) {
					do_write_centre_score( ctr_b, ctr_a, score );
				}

				/// \todo Make these max/min calculations less ugly
				///       (perhaps with C++11 lambdas and nested Boost Range 1.43.0 max/min algorithm functions?)
				min_score = min(min_score, score);
				max_score = max(max_score, score);
			}
		}
		do_plot_scores( scores, min_score, max_score );
	}

}

/// \brief TODOCUMENT
void matrix_plotter::plot_return_path_matrix(const return_path_matrix &arg_return_path_matrix ///< TODOCUMENT
                                             ) {
	check_lengths_match( arg_return_path_matrix.get_length_a(), arg_return_path_matrix.get_length_b() );
	check_window_width_matches( arg_return_path_matrix.get_window_width() );

	const size_t length_a     = get_length_a();
	const size_t length_b     = get_length_b();
	const size_t window_width = get_window_width();

	for (size_t ctr_b = 0; ctr_b < length_b; ++ctr_b) {

		/// \todo Tidy up this window stuff:
		///        * remove the need to worry about offset 1 here
		///        * provide a function that returns the start/stop as a pair
		///        * move towards abstracting iterating over windowed matrices into an iterator of windowed_matrix<>
		const size_t window_start_a = get_window_start_a_for_b__offset_1(
			length_a,
			length_b,
			window_width,
			ctr_b + 1
		) - 1;
		const size_t window_stop_a = get_window_stop_a_for_b__offset_1(
			length_a,
			length_b,
			window_width,
			ctr_b + 1
		) - 1;
		for (size_t ctr_a = window_start_a; ctr_a <= window_stop_a; ++ctr_a) {
			const path_step      next_step  = arg_return_path_matrix.get_path_step_towards_end_at_point( ctr_a, ctr_b );
			const size_size_pair next_point = indices_of_point_after_path_step(next_step, ctr_a, ctr_b);

			// (pairs flipped due to switching from matrix order (row index then column index) to graph order (x then y))
			do_plot_minor_corner_arrow( ctr_b, ctr_a, next_point.second, next_point.first );
		}
	}
}

/// \brief TODOCUMENT
void matrix_plotter::plot_accumulated_scores(const score_accumulation_matrix &arg_score_accumulation_matrix ///< TODOCUMENT
                                             ) {
	check_lengths_match( arg_score_accumulation_matrix.get_length_a(), arg_score_accumulation_matrix.get_length_b() );
	check_window_width_matches( arg_score_accumulation_matrix.get_window_width() );

	const size_t length_a     = get_length_a();
	const size_t length_b     = get_length_b();
	const size_t window_width = get_window_width();

	for (size_t ctr_b = 0; ctr_b < length_b; ++ctr_b) {

		/// \todo Tidy up this window stuff:
		///        * remove the need to worry about offset 1 here
		///        * provide a function that returns the start/stop as a pair
		///        * move towards abstracting iterating over windowed matrices into an iterator of windowed_matrix<>
		const size_t window_start_a = get_window_start_a_for_b__offset_1(
			length_a,
			length_b,
			window_width,
			ctr_b + 1
		) - 1;
		const size_t window_stop_a = get_window_stop_a_for_b__offset_1(
			length_a,
			length_b,
			window_width,
			ctr_b + 1
		) - 1;
		for (size_t ctr_a = window_start_a; ctr_a <= window_stop_a; ++ctr_a) {
			const score_type score = arg_score_accumulation_matrix.get_score_towards_end_at_point( ctr_a, ctr_b );
			// (pair flipped due to switching from matrix order (row index then column index) to graph order (x then y))
			do_write_corner_score( ctr_b, ctr_a, numeric_cast<double>(score) );
		}
	}
}

/// \brief TODOCUMENT
void matrix_plotter::finish(const path &arg_output_stem) const {
	return do_finish(arg_output_stem);
}
