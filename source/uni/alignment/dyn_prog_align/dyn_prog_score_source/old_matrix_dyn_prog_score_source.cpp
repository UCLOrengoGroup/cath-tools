/// \file
/// \brief The old_matrix_dyn_prog_score_source class definitions

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

#include "old_matrix_dyn_prog_score_source.hpp"

#include "ssap/windowed_matrix.hpp"
//#include "common/exception/not_implemented_exception.hpp" // ***** TEMPORARY *****

using namespace cath;
using namespace cath::align;
using namespace std;

using boost::numeric_cast;

/// \brief Return the number of elements in the first entry to by aligned with dynamic-programming
size_t old_matrix_dyn_prog_score_source::do_get_length_a() const {
	return length_a;
}

/// \brief Return the number of elements in the second entry to by aligned with dynamic-programming
size_t old_matrix_dyn_prog_score_source::do_get_length_b() const {
	return length_b;
}

/// \brief TODOCUMENT
score_type old_matrix_dyn_prog_score_source::do_get_score(const size_t &arg_index_a, ///< The index of the element of interest in the first  sequence
                                                          const size_t &arg_index_b  ///< The index of the element of interest in the second sequence
                                                          ) const {
	const int a_matrix_idx = get_window_matrix_a_index__offset_1(length_a, length_b, window, arg_index_a + 1, arg_index_b + 1 );
	return matrix.get().get( arg_index_b + 1, numeric_cast<size_t>( a_matrix_idx ) );
}

/// \brief Ctor for old_matrix_dyn_prog_score_source
old_matrix_dyn_prog_score_source::old_matrix_dyn_prog_score_source(const score_vec_of_vec &arg_matrix,   ///< TODOCUMENT
                                                                   const size_t           &arg_length_a, ///< TODOCUMENT
                                                                   const size_t           &arg_length_b, ///< TODOCUMENT
                                                                   const size_t           &arg_window    ///< TODOCUMENT
                                                                   ) : matrix   ( arg_matrix   ),
                                                                       length_a ( arg_length_a ),
                                                                       length_b ( arg_length_b ),
                                                                       window   ( arg_window   ) {
}

