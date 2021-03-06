/// \file
/// \brief The score_accumulation_matrix class definitions

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

#include "score_accumulation_matrix.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/alignment/dyn_prog_align/detail/return_path_matrix.hpp"
#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/ssap/windowed_matrix.hpp"

using namespace ::cath;
using namespace ::cath::align::detail;
using namespace ::cath::align::gap;
using namespace ::cath::common;
using namespace ::std;

using ::boost::lexical_cast;
using ::boost::numeric_cast;

/// \brief TODOCUMENT
void score_accumulation_matrix::check_length(const size_type &prm_length ///< TODOCUMENT
                                             ) {
	if (prm_length == 0) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Length for score_accumulation_matrix cannot be 0"));
	}
}

/// \brief TODOCUMENT
void score_accumulation_matrix::check_index_against_length(const size_type &prm_index,           ///< TODOCUMENT
                                                           const size_type &prm_length,          ///< TODOCUMENT
                                                           const bool      &prm_permit_one_extra ///< TODOCUMENT
                                                           ) {
	if ( prm_index >= prm_length) {
		if ( prm_index > prm_length || !prm_permit_one_extra ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"score_accumulation index "
				+ lexical_cast<string>(prm_index)
				+ " is too long because corresponding length is "
				+ lexical_cast<string>(prm_length)
			));
		}
	}
}

/// \brief TODOCUMENT
void score_accumulation_matrix::check_indices(const size_type &prm_index_a,         ///< TODOCUMENT
                                              const size_type &prm_index_b,         ///< TODOCUMENT
                                              const bool      &prm_permit_one_extra ///< TODOCUMENT
                                              ) const {
	check_index_against_length(prm_index_a, get_length_a(), prm_permit_one_extra);
	check_index_against_length(prm_index_b, get_length_b(), prm_permit_one_extra);
}

/// \brief TODOCUMENT
void score_accumulation_matrix::initialise(const size_type &prm_length_a,    ///< TODOCUMENT
                                           const size_type &prm_length_b,    ///< TODOCUMENT
                                           const size_type &prm_window_width ///< TODOCUMENT
                                           ) {
	check_length( prm_length_a );
	check_length( prm_length_b );
	check_lengths_and_window_size_are_valid(
		prm_length_a,
		prm_length_b,
		prm_window_width
	);
	scores.assign( prm_length_a, prm_length_b, 0);
	window_width = prm_window_width;
}

/// \brief Ctor for score_accumulation_matrix
score_accumulation_matrix::score_accumulation_matrix(const size_type &prm_length_a,    ///< TODOCUMENT
                                                     const size_type &prm_length_b,    ///< TODOCUMENT
                                                     const size_type &prm_window_width ///< TODOCUMENT
                                                     ) {
	initialise( prm_length_a, prm_length_b, prm_window_width );
}

/// \brief TODOCUMENT
void score_accumulation_matrix::reset(const size_type &prm_length_a,    ///< TODOCUMENT
                                      const size_type &prm_length_b,    ///< TODOCUMENT
                                      const size_type &prm_window_width ///< TODOCUMENT
                                      ) {
	initialise( prm_length_a, prm_length_b, prm_window_width );
}

/// \brief TODOCUMENT
auto score_accumulation_matrix::get_length_a() const -> size_type {
	return scores.get_length_a();
}

/// \brief TODOCUMENT
auto score_accumulation_matrix::get_length_b() const -> size_type {
	return scores.get_length_b();
}

/// \brief TODOCUMENT
auto score_accumulation_matrix::get_window_width() const -> size_type {
	return window_width;
}

/// \brief TODOCUMENT
void score_accumulation_matrix::set_score_towards_end_at_point(const size_type  &prm_index_a, ///< TODOCUMENT
                                                               const size_type  &prm_index_b, ///< TODOCUMENT
                                                               const score_type &prm_score    ///< TODOCUMENT
                                                               ) {
	check_indices( prm_index_a, prm_index_b );
	scores.set( prm_index_a, prm_index_b, prm_score );
}

/// \brief TODOCUMENT
score_type score_accumulation_matrix::get_score_towards_end_at_point(const size_type &prm_index_a, ///< TODOCUMENT
                                                                     const size_type &prm_index_b  ///< TODOCUMENT
                                                                     ) const {
	check_indices(prm_index_a, prm_index_b, true);
	const size_type length_a = get_length_a();
	const size_type length_b = get_length_b();

	// TODOCUMENT
	if ( prm_index_a < length_a && prm_index_b < length_b ) {
		return scores.get( prm_index_a, prm_index_b );
	}
	return 0;
}

/// \brief TODOCUMENT
///
/// \relates score_accumulation_matrix
score_accumulation_matrix cath::align::detail::make_uninitialised_score_accumulation_matrix() {
	return score_accumulation_matrix(1, 1, 1);
}

/// \brief TODOCUMENT
///
/// \relates score_accumulation_matrix
score_type cath::align::detail::get_score_towards_end_from_point_plus_path_step(const score_accumulation_matrix            &prm_score_accumulation_matrix, ///< TODOCUMENT
                                                                                const score_accumulation_matrix::size_type &prm_index_a,                   ///< TODOCUMENT
                                                                                const score_accumulation_matrix::size_type &prm_index_b,                   ///< TODOCUMENT
                                                                                const path_step                            &prm_path_step                  ///< TODOCUMENT
                                                                                ) {
	const size_size_pair path_step_indices = indices_of_point_after_path_step(prm_path_step, prm_index_a, prm_index_b);
	return prm_score_accumulation_matrix.get_score_towards_end_at_point( path_step_indices.first, path_step_indices.second );
}

/// \brief TODOCUMENT
///
/// \relates score_accumulation_matrix
///
/// \relates return_path_matrix
path_step_score_map cath::align::detail::get_total_scores_of_path_steps_from_point(const score_accumulation_matrix &prm_score_accumulation_matrix, ///< TODOCUMENT
                                                                                   const return_path_matrix        &prm_return_path_matrix,        ///< TODOCUMENT
                                                                                   const gap_penalty               &prm_gap_penalty,               ///< TODOCUMENT
                                                                                   const dyn_prog_score_source     &prm_scorer,                    ///< TODOCUMENT
                                                                                   const size_t                    &prm_index_a,                   ///< TODOCUMENT
                                                                                   const size_t                    &prm_index_b                    ///< TODOCUMENT
                                                                                   ) {
	path_step_score_map score_of_path_step;
	for (const path_step &the_path_step : ALL_PATH_STEPS) {
		const score_type step_gap_penalty = get_gap_penalty_for_path_step_from_point(
			prm_return_path_matrix,
			prm_gap_penalty,
			the_path_step,
			prm_index_a,
			prm_index_b
		);
		const score_type beyond_step_score = get_score_towards_end_from_point_plus_path_step(
			prm_score_accumulation_matrix,
			prm_index_a,
			prm_index_b,
			the_path_step
		);

		const score_type step_score = (the_path_step == path_step::ALIGN_PAIR) ? prm_scorer.get_score(prm_index_a, prm_index_b)
		                                                            : numeric_cast<score_type>(0.0);

		const score_type total_score = beyond_step_score + step_score - step_gap_penalty;
		score_of_path_step[the_path_step] = total_score;

//		cerr << prm_index_a;
//		cerr << ", ";
//		cerr << prm_index_b;
//		cerr << " ";
//		cerr << the_path_step;
//		cerr << " : ";
//		cerr << total_score;
//		cerr << " = ";
//		cerr << beyond_step_score;
//		cerr << " + ";
//		cerr << step_score;
//		cerr << " - ";
//		cerr << step_gap_penalty;
//		cerr << endl;
	}
	return score_of_path_step;
}
