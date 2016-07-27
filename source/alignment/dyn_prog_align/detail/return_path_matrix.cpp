/// \file
/// \brief The return_path_matrix class definitions

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

#include <boost/algorithm/string/join.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include "alignment/alignment.h"
#include "alignment/dyn_prog_align/detail/return_path_matrix.h"
#include "alignment/gap/gap_penalty.h"
//#include "common/algorithm/copy_build.h"
#include "common/algorithm/transform_build.h"
#include "exception/invalid_argument_exception.h"
#include "ssap/windowed_matrix.h"

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::align::gap;
using namespace cath::common;
using namespace std;

using boost::adaptors::map_values;
using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::numeric_cast;
using boost::range::max_element;

/// \brief TODOCUMENT
void return_path_matrix::check_length(const size_type &arg_length ///< TODOCUMENT
                                      ) {
	if (arg_length == 0) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Length for return_path_matrix cannot be 0"));
	}
}

/// \brief TODOCUMENT
///
/// Exception guarantee: basic (because one of the assigns in the middle of a list may throw if it can't allocate enough memory)
void return_path_matrix::initialise(const size_type &arg_length_a,    ///< TODOCUMENT
                                    const size_type &arg_length_b,    ///< TODOCUMENT
                                    const size_type &arg_window_width ///< TODOCUMENT
                                    ) {
	check_length( arg_length_a );
	check_length( arg_length_b );
	check_lengths_and_window_size_are_valid(
		arg_length_a,
		arg_length_b,
		arg_window_width
	);

	return_path.resize( arg_length_a );
	for (path_step_vec &return_path_line : return_path) {
		return_path_line.assign( arg_length_b, path_step::ALIGN_PAIR );
	}
	window_width = arg_window_width;
}

/// \brief Ctor for return_path_matrix
return_path_matrix::return_path_matrix(const size_type &arg_length_a,    ///< TODOCUMENT
                                       const size_type &arg_length_b,    ///< TODOCUMENT
                                       const size_type &arg_window_width ///< TODOCUMENT
                                       ) {
	initialise( arg_length_a, arg_length_b, arg_window_width );
}

/// \brief TODOCUMENT
void return_path_matrix::reset(const size_type &arg_length_a,    ///< TODOCUMENT
                               const size_type &arg_length_b,    ///< TODOCUMENT
                               const size_type &arg_window_width ///< TODOCUMENT
                               ) {
	initialise( arg_length_a, arg_length_b, arg_window_width );
}

/// \brief TODOCUMENT
return_path_matrix::size_type return_path_matrix::get_length_a() const {
	return return_path.size();
}


/// \brief TODOCUMENT
return_path_matrix::size_type return_path_matrix::get_length_b() const {
	return return_path.front().size();
}

/// \brief TODOCUMENT
return_path_matrix::size_type return_path_matrix::get_window_width() const {
	return window_width;
}

/// \brief TODOCUMENT
void return_path_matrix::set_path_step_towards_end_at_point(const size_type &arg_index_of_next_element_in_a, ///< TODOCUMENT
                                                            const size_type &arg_index_of_next_element_in_b, ///< TODOCUMENT
                                                            const path_step &arg_path_step                   ///< TODOCUMENT
                                                            ) {
	check_indices_are_within_window(
		get_length_a(),
		get_length_b(),
		window_width,
		arg_index_of_next_element_in_a,
		arg_index_of_next_element_in_b
	);
	return_path[arg_index_of_next_element_in_a][arg_index_of_next_element_in_b] = arg_path_step;
}

/// \brief TODOCUMENT
path_step return_path_matrix::get_path_step_towards_end_at_point(const size_type &arg_index_of_next_element_in_a, ///< TODOCUMENT
                                                                 const size_type &arg_index_of_next_element_in_b  ///< TODOCUMENT
                                                                 ) const {
	const size_type length_a = get_length_a();
	const size_type length_b = get_length_b();
	if ( arg_index_of_next_element_in_a <  length_a && arg_index_of_next_element_in_b == length_b ) {
		return path_step::INSERT_INTO_FIRST;
	}
	if ( arg_index_of_next_element_in_a == length_a && arg_index_of_next_element_in_b <  length_b ) {
		return path_step::INSERT_INTO_SECOND;
	}

	check_indices_are_within_window(
		length_a,
		length_b,
		window_width,
		arg_index_of_next_element_in_a,
		arg_index_of_next_element_in_b
	);
	return return_path[arg_index_of_next_element_in_a][arg_index_of_next_element_in_b];
}

/// \brief TODOCUMENT
///
/// \relates return_path_matrix
return_path_matrix cath::align::detail::make_uninitialised_return_path_matrix() {
	return return_path_matrix(1, 1, 1);
}

/// \brief Build an alignment from a return path matrix by tracing the path from the start (top-left)
///
/// Exception guarantee: strong
///
/// \relates return_path_matrix
///
/// \relates alignment
alignment cath::align::detail::make_alignment(const return_path_matrix &arg_return_path_matrix ///< TODOCUMENT
                                              ) {
	// Grab the lengths of the two entities represented by the return_path_matrix
	const return_path_matrix::size_type length_a = arg_return_path_matrix.get_length_a();
	const return_path_matrix::size_type length_b = arg_return_path_matrix.get_length_b();

	// Create a new alignment and reserve enough memory to at least store the length of the longer
	alignment new_alignment( alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT );
	new_alignment.reserve( max( length_a, length_b ) );

	// Start from (0, 0) (the top-left) and keep iterating until hitting the end of both (the bottom-right)
	size_size_pair position(0, 0);
	while (position.first < length_a || position.second < length_b) {
		// Grab the next step from this position
		const path_step next_step = arg_return_path_matrix.get_path_step_towards_end_at_point(
			position.first,
			position.second
		);

		// Append this step to the alignment
		append_path_step_to_pair_alignment_from_point(new_alignment, next_step, position);

		// Update the position
		position = indices_of_point_after_path_step(next_step, position.first, position.second);
	}

//	// If necessary, finish by appending the entries corresponding to moving along the
//	// bottom or right edge of the matrix to get to the bottom-right corner
//	// (the return_path_matrix doesn't store these path_steps because they're just implied)
//	while (position.first < length_a) {
//		append_path_step_to_pair_alignment_from_point(new_alignment, path_step::INSERT_INTO_FIRST, position);
//		++position.first;
//	}
//	while (position.second < length_b) {
//		append_path_step_to_pair_alignment_from_point(new_alignment, path_step::INSERT_INTO_SECOND, position);
//		++position.second;
//	}

	// Return the resulting alignment
	return new_alignment;
}

/// \brief TODOCUMENT
///
/// \relates return_path_matrix
path_step cath::align::detail::get_path_dirn_towards_end_from_point_plus_path_step(const return_path_matrix            &arg_return_path_matrix, ///< TODOCUMENT
                                                                                   const return_path_matrix::size_type &arg_index_a,            ///< TODOCUMENT
                                                                                   const return_path_matrix::size_type &arg_index_b,            ///< TODOCUMENT
                                                                                   const path_step                     &arg_path_step           ///< TODOCUMENT
                                                                                   ) {
	const size_size_pair path_step_indices = indices_of_point_after_path_step(arg_path_step, arg_index_a, arg_index_b);
	return arg_return_path_matrix.get_path_step_towards_end_at_point( path_step_indices.first, path_step_indices.second );
}

/// \brief TODOCUMENT
///
/// \relates return_path_matrix
score_type cath::align::detail::get_gap_penalty_for_path_step_from_point(const return_path_matrix            &arg_return_path_matrix, ///< TODOCUMENT
                                                                         const gap_penalty                   &arg_gap_penalty,        ///< TODOCUMENT
                                                                         const path_step                     &arg_path_step,          ///< TODOCUMENT
                                                                         const return_path_matrix::size_type &arg_index_a,            ///< TODOCUMENT
                                                                         const return_path_matrix::size_type &arg_index_b             ///< TODOCUMENT
                                                                         ) {
	const score_type zero_score = numeric_cast<score_type>(0.0);
	// If this is an aligned pair, then there is no gap here so return zero
	if ( arg_path_step == path_step::ALIGN_PAIR                             ) {
		return zero_score;
	}
	// Or if this is an insert into the first and this is still at the start of the second, then this isn't a gap so return 0
	if ( arg_path_step == path_step::INSERT_INTO_FIRST  && arg_index_b == 0 ) {
		return zero_score;
	}
	// Or if this is an insert into the second and this is still at the start of the first, then this isn't a gap so return 0
	if ( arg_path_step == path_step::INSERT_INTO_SECOND && arg_index_a == 0 ) {
		return zero_score;
	}
	const path_step path_step_after_path_step = get_path_dirn_towards_end_from_point_plus_path_step(
		arg_return_path_matrix,
		arg_index_a,
		arg_index_b,
		arg_path_step
	);
	return (path_step_after_path_step == arg_path_step) ? arg_gap_penalty.get_extend_gap_penalty()
	                                                    : arg_gap_penalty.get_open_gap_penalty();
}

/// \brief TODOCUMENT
///
/// \relates return_path_matrix
size_size_pair cath::align::detail::get_b_window_start_and_stop_for_a_index(const return_path_matrix            &arg_return_path_matrix, ///< TODOCUMENT
                                                                            const return_path_matrix::size_type &arg_index_a             ///< TODOCUMENT
                                                                            ) {
	const return_path_matrix::size_type length_a       = arg_return_path_matrix.get_length_a();
	const return_path_matrix::size_type length_b       = arg_return_path_matrix.get_length_b();
	const return_path_matrix::size_type window_width   = arg_return_path_matrix.get_window_width();
	const size_t window_start_b = get_window_start_a_for_b__offset_1(
		length_b,
		length_a,
		window_width,
		arg_index_a + 1
	) - 1;
	const size_t window_stop_b  = get_window_stop_a_for_b__offset_1(
		length_b,
		length_a,
		window_width,
		arg_index_a + 1
	) - 1;
	return make_pair(window_start_b, window_stop_b);
}

/// \brief TODOCUMENT
///
/// \relates return_path_matrix
ostream & cath::align::detail::operator<<(ostream                  &arg_os,                ///< TODOCUMENT
                                          const return_path_matrix &arg_return_path_matrix ///< TODOCUMENT
                                          ) {
	using path_size_type = return_path_matrix::size_type;

	const path_size_type length_a = arg_return_path_matrix.get_length_a();
	const path_size_type length_b = arg_return_path_matrix.get_length_b();
	for (path_size_type ctr_a = 0; ctr_a < length_a; ++ctr_a) {
		const size_size_pair b_window_start_and_stop     = get_b_window_start_and_stop_for_a_index(
			arg_return_path_matrix,
			ctr_a
		);
		const size_size_pair window_of_prev_or_first_row = get_b_window_start_and_stop_for_a_index(
			arg_return_path_matrix,
			max(numeric_cast<path_size_type>(1), ctr_a) - 1
		);
		str_str_str_str_tpl lines;
		for (path_size_type ctr_b = 0; ctr_b < length_b; ++ctr_b) {
			const bool within_window = (ctr_b >= b_window_start_and_stop.first && ctr_b <= b_window_start_and_stop.second);
			if ( !within_window ) {
				get<0>( lines ) += ".      ";
				get<1>( lines ) += "       ";
				get<2>( lines ) += "       ";
				get<3>( lines ) += "       ";
			}
			else {
				const path_step the_path_step = arg_return_path_matrix.get_path_step_towards_end_at_point(ctr_a, ctr_b);
				switch (the_path_step) {
					case ( path_step::ALIGN_PAIR         ) : {
						get<0>( lines ) += "@_     ";
						get<1>( lines ) += "  \\__  ";
						get<2>( lines ) += "     \\ ";
						get<3>( lines ) += "      4";
						break;
					}
					case ( path_step::INSERT_INTO_FIRST  ) : {
						get<0>( lines ) += "@      ";
						get<1>( lines ) += "|      ";
						get<2>( lines ) += "|      ";
						get<3>( lines ) += "V      ";
						break;
					}
					case ( path_step::INSERT_INTO_SECOND ) : {
						get<0>( lines ) += "@----->";
						get<1>( lines ) += "       ";
						get<2>( lines ) += "       ";
						get<3>( lines ) += "       ";
						break;
					}
				}
			}
		}

		if ( window_of_prev_or_first_row.second >= length_b - 1 ) {
			get<0>( lines ) += "@";
			get<1>( lines ) += "|";
			get<2>( lines ) += "|";
			get<3>( lines ) += "V";
		}
		else {
			get<0>( lines ) += ".";
			get<1>( lines ) += " ";
			get<2>( lines ) += " ";
			get<3>( lines ) += " ";
		}
		arg_os << get<0>( lines ) << "\n";
		arg_os << get<1>( lines ) << "\n";
		arg_os << get<2>( lines ) << "\n";
		arg_os << get<3>( lines ) << "\n";
	}
	const size_size_pair last_window_start_and_stop = get_b_window_start_and_stop_for_a_index(
		arg_return_path_matrix,
		length_a - 1
	);
	const size_t last_window_width = (last_window_start_and_stop.second - last_window_start_and_stop.first);
	arg_os << join(str_vec(length_b - last_window_width, ".      "), "");
	arg_os << join(str_vec(           last_window_width, "@----->"), "");
	arg_os << "@\n";
	return arg_os;
}

/// \brief TODOCUMENT
///
/// \relates path_step
score_type cath::align::detail::max_path_step_score(const path_step_score_map &arg_score_of_path_step ///< TODOCUMENT
                                                    ) {
	return *max_element(
		arg_score_of_path_step | map_values
	);
}
