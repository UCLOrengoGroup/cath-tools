/// \file
/// \brief The superpose_orderer class definitions

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

#include "superpose_orderer.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "common/algorithm/sort_copy.h"
#include "common/c++14/cbegin_cend.h"
#include "exception/invalid_argument_exception.h"
#include "exception/not_implemented_exception.h"
#include "exception/runtime_error_exception.h"
#include "superposition/detail/spanning_tree_greater.h"

#include <algorithm>

using namespace cath;
using namespace cath::common;
using namespace cath::sup;
using namespace cath::sup::detail;
using namespace std;

using boost::adjacency_list;
//using boost::algorithm::join;
//using boost::algorithm::replace_all_copy;
using boost::edge_weight_t;
using boost::graph_traits;
using boost::kruskal_minimum_spanning_tree;
using boost::lexical_cast;
using boost::no_property;
using boost::property;
using boost::range::sort;
using boost::undirectedS;
using boost::vecS;

/// \brief TODOCUMENT
superpose_orderer::size_type superpose_orderer::half_matrix_index_of_indices(const size_type &arg_row_index, ///< TODOCUMENT
                                                                             const size_type &arg_col_index  ///< TODOCUMENT
                                                                             ) {
	// Sanity check that the row index is greater than the column index
	if (arg_row_index <= arg_col_index) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"superpose_orderer: row index "  + lexical_cast<string>( arg_row_index )
			+ " should exceed column index " + lexical_cast<string>( arg_col_index )
		));
	}
	return ( ( arg_row_index * ( arg_row_index -1 ) / 2 ) + arg_col_index );
}

/// \brief TODOCUMENT
void superpose_orderer::check_index_pair(const size_type &arg_row_index, ///< TODOCUMENT
                                         const size_type &arg_col_index  ///< TODOCUMENT
                                         ) const {
	// Check that both indices are within the number of items
	if (arg_row_index >= get_num_items() || arg_col_index >= get_num_items()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("superpose_orderer: index is out of range"));
	}

	// Check that the relationship between the row_index and col_index is valid
	// (and double check that they point to a slot with in scores
	const size_type half_matrix_index = half_matrix_index_of_indices(arg_row_index, arg_col_index);
	if (half_matrix_index >= scores.size()) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("superpose_orderer: row and column indices appeared valid but they overrun the allocated memory"));
	}
}

/// \brief Ctor for superpose_orderer
superpose_orderer::superpose_orderer(const size_type &arg_num_items ///< TODOCUMENT
                                     ) : num_items(arg_num_items) {
	if (get_num_items() < 1) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("superpose_orderer must have some strictly positive number of items"));
	}
	const size_type half_matrix_size = half_matrix_index_of_indices(get_num_items(), 0);
	is_scored.assign( half_matrix_size, false );
	scores.assign(    half_matrix_size, 0.0   );
}

/// \brief TODOCUMENT
superpose_orderer::size_type superpose_orderer::get_num_items() const {
	return num_items;
}

/// \brief TODOCUMENT
bool superpose_orderer::has_score(const size_type &arg_row_index, ///< TODOCUMENT
                                  const size_type &arg_col_index  ///< TODOCUMENT
                                  ) const {
	check_index_pair(arg_row_index, arg_col_index);
	const size_type half_matrix_index = half_matrix_index_of_indices(arg_row_index, arg_col_index);
	return is_scored[half_matrix_index];
}

/// \brief TODOCUMENT
double superpose_orderer::get_score(const size_type &arg_row_index, ///< TODOCUMENT
                                    const size_type &arg_col_index  ///< TODOCUMENT
                                    ) const {
	if (!has_score(arg_row_index, arg_col_index)) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("superpose_orderer::get_score() requested to get score from a location for which no score (yet) exists"));
	}
	const size_type half_matrix_index = half_matrix_index_of_indices(arg_row_index, arg_col_index);
	return scores[half_matrix_index];
}

/// \brief TODOCUMENT
void superpose_orderer::set_score(const size_type &arg_row_index, ///< TODOCUMENT
                                  const size_type &arg_col_index, ///< TODOCUMENT
                                  const double    &arg_score      ///< TODOCUMENT
                                  ) {
	check_index_pair(arg_row_index, arg_col_index);
	if (!boost::math::isfinite(arg_score)) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("superpose_orderer::set_score() requires a score that is a normal number"));
	}
	const size_type half_matrix_index = half_matrix_index_of_indices(arg_row_index, arg_col_index);
	is_scored[ half_matrix_index ] = true;
	scores[    half_matrix_index ] = arg_score;
}

/// \brief Find a maximum spanning tree of the specified superpose_orderer
///
/// \relates superpose_orderer
///
/// \returns A vector of pair<size_t, size_t>, sorted by the standard pair ordering
///          (ie ascending on the first value and then on the second)
///
/// This is implemented by using Boost Graph library's kruskal_minimum_spanning_tree()
/// so much of this subroutine involves preparing data for that call and then extracting the results
size_size_pair_vec cath::sup::get_spanning_tree(const superpose_orderer &arg_superpose_orderer ///< The superpose_orderer for which to find a maximum spanning tree
                                                ) {
	// Prepare some type aliases that are useful for this
	using Graph = adjacency_list < vecS, vecS, undirectedS, no_property, property <edge_weight_t, double> >;
	using edge_desc = graph_traits < Graph >::edge_descriptor;

	// Construct edges and weights vectors and reserve enough memory for the half matrix
	const superpose_orderer::size_type num_items = arg_superpose_orderer.get_num_items();
	size_size_pair_vec edges;
	edges.reserve((num_items * (num_items-1)) / 2);
	doub_vec weights;
	weights.reserve((num_items * (num_items-1)) / 2);

	// Populate the edges and weights vector with the data from arg_superpose_orderer.
	// Use the negative of the score so that kruskal_minimum_spanning_tree() finds the maximum spanning tree.
	for (superpose_orderer::size_type row_ctr = 0; row_ctr < num_items; ++row_ctr) {
		for (superpose_orderer::size_type col_ctr = 0; col_ctr < num_items; ++col_ctr) {
			if (col_ctr < row_ctr) {
				if (arg_superpose_orderer.has_score(row_ctr, col_ctr)) {
					edges.push_back(make_pair(row_ctr, col_ctr));
					weights.push_back( - arg_superpose_orderer.get_score(row_ctr, col_ctr));
				}
			}
		}
	}

	// Construct a graph from the edges and weights
	const Graph my_graph(
		common::cbegin( edges   ),
		common::cend  ( edges   ),
		common::cbegin( weights ),
		num_items
	);

	// Call kruskal_minimum_spanning_tree() to construct the spanning tree
	vector<edge_desc> spanning_tree;
	kruskal_minimum_spanning_tree(my_graph, back_inserter(spanning_tree));

	// Convert the results into vector of size_size_pairs
	// The target is given first since this will be the lower value
	size_size_pair_vec pairs;
	pairs.reserve(spanning_tree.size());
	for (const edge_desc &spanning_tree_edge : spanning_tree) {
		pairs.push_back(make_pair(target(spanning_tree_edge, my_graph), source(spanning_tree_edge, my_graph)));
	}

	// If the number of edges is not one less than the number of items, then it was not
	// possible to find a single tree to span all items so throw an exception
	if (pairs.size() + 1 != num_items) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("get_order_of_pairs_to_superpose() : Unable to find tree to span all items"));
	}

	// Return a sorted copy of the pairs
	return sort_copy( pairs );
}

/// \brief TODOCUMENT
///
/// \relates superpose_orderer
size_size_pair_vec cath::sup::get_spanning_tree_ordered_by_desc_score(const superpose_orderer &arg_superpose_orderer ///< TODOCUMENT
                                                                      ) {
	const size_size_pair_vec orig_spanning_tree = get_spanning_tree( arg_superpose_orderer );
	const size_size_pair_vec new_spanning_tree  = order_spanning_tree_by_desc_score( arg_superpose_orderer, orig_spanning_tree );
	return new_spanning_tree;
}

/// \brief TODOCUMENT
///
/// \relates superpose_orderer
size_size_pair_vec cath::sup::order_spanning_tree_by_desc_score(const superpose_orderer  &arg_superpose_orderer, ///< TODOCUMENT
                                                                const size_size_pair_vec &arg_spanning_tree      ///< TODOCUMENT
                                                                ) {
	size_size_pair_vec new_spanning_tree = arg_spanning_tree;
	sort(
		new_spanning_tree,
		spanning_tree_greater( arg_superpose_orderer )
	);
	return new_spanning_tree;
}

/// \brief TODOCUMENT
///
/// \relates superpose_orderer
superpose_orderer cath::sup::make_superpose_orderer(const size_size_pair_doub_map &arg_score_by_index_pair ///< TODOCUMENT
                                                    ) {
	// Find the largest index in arg_score_by_index_pair
	size_t largest_index = 0;
	for (const size_size_pair_doub_map_value &index_pair_and_score : arg_score_by_index_pair) {
		largest_index = max( largest_index, index_pair_and_score.first.first  );
		largest_index = max( largest_index, index_pair_and_score.first.second );
	}

	// Load the scores into the new orderer
	superpose_orderer blue_monday( largest_index + 1 );
	for (const size_size_pair_doub_map_value &index_pair_and_score : arg_score_by_index_pair) {
		const size_t &index_1    = index_pair_and_score.first.first;
		const size_t &index_2    = index_pair_and_score.first.second;
		const double &ssap_score = index_pair_and_score.second;
		blue_monday.set_score( max(index_1, index_2), min(index_1, index_2), ssap_score );
	}
	return blue_monday;
}

/// \brief TODOCUMENT
///
/// \relates superpose_orderer
size_size_pair_vec cath::sup::get_spanning_tree_ordered_by_desc_score(const size_size_pair_doub_map &arg_score_by_index_pair ///< TODOCUMENT
                                                                      ) {
	const superpose_orderer the_sup_orderer = make_superpose_orderer( arg_score_by_index_pair );
	return get_spanning_tree_ordered_by_desc_score( the_sup_orderer );
}

