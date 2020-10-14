/// \file
/// \brief The new_matrix_dyn_prog_score_source class definitions

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

#include "new_matrix_dyn_prog_score_source.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include "ssap/windowed_matrix.hpp"

using namespace cath;
using namespace cath::align;
using namespace std;
using boost::numeric_cast;

/// \brief Return the number of elements in the first entry to by aligned with dynamic-programming
size_t new_matrix_dyn_prog_score_source::do_get_length_a() const {
	return length_a;
}

/// \brief Return the number of elements in the second entry to by aligned with dynamic-programming
size_t new_matrix_dyn_prog_score_source::do_get_length_b() const {
	return length_b;
}

/// \brief TODOCUMENT
score_type new_matrix_dyn_prog_score_source::do_get_score(const size_t &prm_index_a, ///< The index of the element of interest in the first  sequence
                                                          const size_t &prm_index_b  ///< The index of the element of interest in the second sequence
                                                          ) const {
	return numeric_cast<score_type>(matrix[ prm_index_a ][ prm_index_b ]);
}

/// \brief Ctor for new_matrix_dyn_prog_score_source
new_matrix_dyn_prog_score_source::new_matrix_dyn_prog_score_source(const float_score_vec_vec &prm_matrix,   ///< TODOCUMENT
                                                                   const size_t              &prm_length_a, ///< TODOCUMENT
                                                                   const size_t              &prm_length_b  ///< TODOCUMENT
                                                                   ) : matrix   ( prm_matrix   ),
                                                                       length_a ( prm_length_a ),
                                                                       length_b ( prm_length_b ) {
}
