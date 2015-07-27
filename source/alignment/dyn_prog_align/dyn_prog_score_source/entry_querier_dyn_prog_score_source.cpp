/// \file
/// \brief The entry_querier_dyn_prog_score_source class definitions

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

#include "entry_querier_dyn_prog_score_source.h"

#include "structure/entry_querier/entry_querier.h"

using namespace cath;
using namespace cath::align;
//using namespace std;

/// \brief Return the number of elements in the first entry to by aligned with dynamic-programming
size_t entry_querier_dyn_prog_score_source::do_get_length_a() const {
	return the_entry_querier.get_length(protein_a);
}

/// \brief Return the number of elements in the second entry to by aligned with dynamic-programming
size_t entry_querier_dyn_prog_score_source::do_get_length_b() const {
	return the_entry_querier.get_length(protein_b);
}

/// \brief TODOCUMENT
score_type entry_querier_dyn_prog_score_source::do_get_score(const size_t &arg_index_a, ///< The index of the element of interest in the first  sequence
                                                             const size_t &arg_index_b  ///< The index of the element of interest in the second sequence
                                                             ) const {
	return the_entry_querier.distance_score__offset_1(
		protein_a,             protein_b,
		view_from_index_a + 1, view_from_index_b + 1,
		arg_index_a       + 1, arg_index_b       + 1
	);
}

/// \brief Ctor for entry_querier_dyn_prog_score_source
entry_querier_dyn_prog_score_source::entry_querier_dyn_prog_score_source(const entry_querier &arg_entry_querier,     ///< TODOCUMENT
                                                                         const protein       &arg_protein_a,         ///< TODOCUMENT
                                                                         const protein       &arg_protein_b,         ///< TODOCUMENT
                                                                         const size_t        &arg_view_from_index_a, ///< TODOCUMENT
                                                                         const size_t        &arg_view_from_index_b  ///< TODOCUMENT
                                                                         ) : the_entry_querier ( arg_entry_querier     ),
                                                                             protein_a         ( arg_protein_a         ),
                                                                             protein_b         ( arg_protein_b         ),
                                                                             view_from_index_a ( arg_view_from_index_a ),
                                                                             view_from_index_b ( arg_view_from_index_b ) {
}
