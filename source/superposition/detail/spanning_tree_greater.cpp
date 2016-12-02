/// \file
/// \brief The spanning_tree_greater class definitions

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

#include "spanning_tree_greater.hpp"

#include "superposition/superpose_orderer.hpp"

using namespace cath::sup::detail;
using namespace std;

/// \brief Ctor for spanning_tree_greater
spanning_tree_greater::spanning_tree_greater(const superpose_orderer &arg_superpose_orderer ///< The superpose_orderer from which to get the scores for pairs of indices
                                             ) : the_superpose_orderer(arg_superpose_orderer) {
}

/// \brief Return whether the superpose_orderer's score for the first pair is greater than that of the second pair
bool spanning_tree_greater::operator()(const size_size_pair &arg_index_a, ///< The first  pair of indices to compare
                                       const size_size_pair &arg_index_b  ///< The second pair of indices to compare
                                       ) const {
	const double score_a = the_superpose_orderer.get_score(
		max( arg_index_a.first, arg_index_a.second ),
		min( arg_index_a.first, arg_index_a.second )
	);
	const double score_b = the_superpose_orderer.get_score(
		max( arg_index_b.first, arg_index_b.second ),
		min( arg_index_b.first, arg_index_b.second )
	);
	return ( score_a > score_b );
}
