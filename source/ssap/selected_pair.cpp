/// \file
/// \brief The selected_pair class definitions

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

#include "selected_pair.h"

using namespace cath;

/// \brief Simple ctor
selected_pair::selected_pair(const size_t     &arg_index_a, ///< The pair's first  index
                             const size_t     &arg_index_b, ///< The pair's second index
                             const score_type &arg_score    ///< The pair's score
                             ) : index_a(arg_index_a),
                                 index_b(arg_index_b),
                                 score  (arg_score) {
}

/// \brief Getter for first index
size_t selected_pair::get_index_a() const {
	return index_a;
}

/// \brief Getter for second index
size_t selected_pair::get_index_b() const {
	return index_b;
}

/// \brief Getter for the score
score_type selected_pair::get_score() const {
	return score;
}

/// \brief Less-than operator that returns the result of a less-than comparison of the selected_pairs' scores
bool cath::operator<(const selected_pair &arg_selected_pair_1, ///< The first  selected_pair to be compared
                     const selected_pair &arg_selected_pair_2  ///< The second selected_pair to be compared
                     ) {
	return (arg_selected_pair_1.get_score() < arg_selected_pair_2.get_score());
}
