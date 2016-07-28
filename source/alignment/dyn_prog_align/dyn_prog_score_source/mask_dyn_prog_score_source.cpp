/// \file
/// \brief The mask_dyn_prog_score_source class definitions

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

#include "mask_dyn_prog_score_source.h"

using namespace cath;
using namespace cath::align;
//using namespace std;

/// \brief Return the number of elements in the first entry to by aligned with dynamic-programming
size_t mask_dyn_prog_score_source::do_get_length_a() const {
	return masked_score_source.get_length_a();
}

/// \brief Return the number of elements in the second entry to by aligned with dynamic-programming
size_t mask_dyn_prog_score_source::do_get_length_b() const {
	return masked_score_source.get_length_b();
}

/// \brief TODOCUMENT
score_type mask_dyn_prog_score_source::do_get_score(const size_t &arg_index_a, ///< The index of the element of interest in the first  sequence
                                                    const size_t &arg_index_b  ///< The index of the element of interest in the second sequence
                                                    ) const {
	// Grab the mask value for this entry from the mask_matrix
	const bool mask_result = mask_matrix.get( arg_index_b + 1, arg_index_a + 1 );

	// If the mask value is true then return the masked_score_source's score, else return 0
	return mask_result ? masked_score_source.get_score(arg_index_a, arg_index_b)
	                   : 0;
}

/// \brief Ctor for mask_dyn_prog_score_source
mask_dyn_prog_score_source::mask_dyn_prog_score_source(const bool_vec_of_vec       &arg_mask_matrix,        ///< TODOCUMENT
                                                       const dyn_prog_score_source &arg_masked_score_source ///< TODOCUMENT
                                                       ) : mask_matrix        ( arg_mask_matrix         ),
                                                           masked_score_source( arg_masked_score_source ) {
}
