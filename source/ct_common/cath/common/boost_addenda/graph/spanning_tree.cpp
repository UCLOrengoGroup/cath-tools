/// \file
/// \brief The spanning_tree class definitions

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

#include "spanning_tree.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <tuple>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/size_t_literal.hpp"

using namespace ::cath;
using namespace ::cath::common;

using ::boost::adaptors::filtered;
using ::boost::adaptors::transformed;
using ::boost::algorithm::join;
using ::boost::irange;
using ::std::filesystem::path;
using ::std::make_pair;
using ::std::make_tuple;
using ::std::max;
using ::std::min;
using ::std::ofstream;
using ::std::string;
using ::std::vector;

/// \brief Get a simple ( 0 <-> 1 <-> 2 <-> ... <-> (prm_num_items -1) spanning tree)
size_size_pair_vec cath::common::make_simple_unweighted_spanning_tree(const size_t &prm_num_items ///< The number of items to span
                                                                      ) {
	// If zero/one items then return empty
	if ( prm_num_items <= 1 ) {
		return {};
	}

	// Otherwise return a spanning tree between adjacent indices
	return transform_build<size_size_pair_vec>(
		irange( 1_z, prm_num_items ),
		[] (const size_t &x) { return make_pair( x - 1, x ); }
	);
}

/// \brief Calculate a max-spanning-tree over the specified edges_and_scores and number of items
size_size_doub_tpl_vec cath::common::calc_max_spanning_tree(const size_size_doub_tpl_vec &prm_edges,    ///< The weighted edges from which a spanning tree should be extracted
                                                            const size_t                 &prm_num_items ///< The number of items to span
                                                            ) {
	return calc_max_spanning_tree(
		prm_edges | transformed( [] (const size_size_doub_tpl &x) { return make_pair( get<0>( x ), get<1>( x ) ); } ),
		prm_edges | transformed( [] (const size_size_doub_tpl &x) { return get<2>( x );                           } ),
		prm_num_items
	);
}

/// \brief Calculate a min-spanning-tree over the specified edges_and_scores and number of items
size_size_doub_tpl_vec cath::common::calc_min_spanning_tree(const size_size_doub_tpl_vec &prm_edges,    ///< The weighted edges from which a spanning tree should be extracted
                                                            const size_t                 &prm_num_items ///< The number of items to span
                                                            ) {
	return calc_min_spanning_tree(
		prm_edges | transformed( [] (const size_size_doub_tpl &x) { return make_pair( get<0>( x ), get<1>( x ) ); } ),
		prm_edges | transformed( [] (const size_size_doub_tpl &x) { return get<2>( x );                           } ),
		prm_num_items
	);
}

/// \brief Get the edges of the specified spanning tree (ie strip off the weights)
size_size_pair_vec cath::common::get_edges_of_spanning_tree(const size_size_doub_tpl_vec &prm_spanning_tree ///< The weighted spanning-tree from which the weights should be stripped
                                                            ) {
	return transform_build<size_size_pair_vec>(
		prm_spanning_tree,
		[] (const size_size_doub_tpl &x) {
			return make_pair( get<0>( x ), get<1>( x ) );
		}
	);
}

/// \brief Return a copy of the specified spanning tree ordered such that:
///          * the first edge is the one with the specified index
///          * all edges after that contain exactly one node contained in an edge before it
///          * within those constraints, edges with higher weight/scores are always added before others
///          * within those constraints, edges with higher highest node index are always added before others
///          * within those constraints, edges with higher lowest  node index are always added before others
size_size_doub_tpl_vec cath::common::order_spanning_tree_from_start(const size_size_doub_tpl_vec &prm_spanning_tree, ///< The spanning tree to process
                                                                    const size_t                 &prm_index          ///< The index of the edge at which to start
                                                                    ) {
	if ( prm_index >= prm_spanning_tree.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot start the spanning tree from an edge with an index >= the number of edges"));
	}

	// Function object for generating the preferability of an edge corresponding to an index
	const auto index_preferability = [&] (const size_t &x) {
		const auto &edge = prm_spanning_tree[ x ];
		return make_tuple(
			get<2>( edge ),
			max( get<0>( edge ), get<1>( edge ) ),
			min( get<0>( edge ), get<1>( edge ) )
		);
	};

	// Build a data structure of the edge indices involved in the spanning tree
	// (excluding prm_index), each stored under the indices of their two nodes
	vector<size_set> edge_indices_by_node = [&] {
		vector<size_set> edge_indcs( prm_spanning_tree.size() + 1 );
		const auto non_starting_edge_indices = indices( prm_spanning_tree.size() )
			| filtered( [&] (const size_t &x) { return x != prm_index; } );
		for (const size_t &non_starting_edge_index : non_starting_edge_indices) {
			const auto &[ node_a, node_b, val ] = prm_spanning_tree[ non_starting_edge_index ];
			edge_indcs[ node_a ].insert( non_starting_edge_index );
			edge_indcs[ node_b ].insert( non_starting_edge_index );
		}
		return edge_indcs;
	} ();

	// Get the first edge
	const size_size_doub_tpl &first_edge = prm_spanning_tree[ prm_index ];

	// Record the nodes included in edges so far
	size_set nodes_so_far = { get<0>( first_edge ), get<1>( first_edge ) };

	// Return a size_size_doub_tpl_vec with one element for each value of a counter
	return transform_build<size_size_doub_tpl_vec>(
		indices( prm_spanning_tree.size() ),
		[&] (const size_t &ctr) {
			// For the first element, just return the first edge
			if ( ctr == 0 ) {
				return first_edge;
			}

			// Otherwise, find the index of the best edge from the nodes visited so far
			size_opt index_of_best;
			for (const size_t &node_so_far : nodes_so_far) {
				for (const size_t &edge_index : edge_indices_by_node[ node_so_far ] ) {
					if ( ! index_of_best || index_preferability( edge_index ) > index_preferability( *index_of_best ) ) {
						index_of_best = edge_index;
					}
				}
			}

			// Update nodes_so_far and edge_indices_by_node with the new best edge
			const auto   &best_edge   = prm_spanning_tree[ *index_of_best ];
			const size_t &best_node_a = get<0>( best_edge );
			const size_t &best_node_b = get<1>( best_edge );

			edge_indices_by_node[ best_node_a ].erase( *index_of_best );
			edge_indices_by_node[ best_node_b ].erase( *index_of_best );

			nodes_so_far.insert( best_node_a );
			nodes_so_far.insert( best_node_b );

			// Return the best edge
			return best_edge;
		}
	);
}

/// \brief Create a graphviz string representing the specified spanning tree
string cath::common::make_graphviz_string_of_spanning_tree(const size_size_doub_tpl_vec &prm_spanning_tree ///< The tree to represent in graphviz format
                                                           ) {
	using ::std::to_string;
	return "digraph example {\n"
		+ join(
			indices( prm_spanning_tree.size() )
				| transformed( [&] (const size_t &x) {
					return
						  to_string( get<0>( prm_spanning_tree[ x ] ) )
						+ " -> "
						+ to_string( get<1>( prm_spanning_tree[ x ] ) )
						+ " [ "
						+ to_string( get<2>( prm_spanning_tree[ x ] ) )
						+ ", "
						+ to_string( x )
						+ "]"
						+ "\n";
				} ),
			""
		)
		+ "}\n";
}

/// \brief Write the specified spanning tree as a graphviz string to the specified file
void cath::common::write_graphviz_string_of_spanning_tree(const size_size_doub_tpl_vec &prm_spanning_tree, ///< The spanning tree to represent in the specified file as a graphviz string
                                                          const path                   &prm_file           ///< The file to write to
                                                          ) {
	ofstream output_stream = open_ofstream( prm_file );
	output_stream << make_graphviz_string_of_spanning_tree( prm_spanning_tree );
	output_stream.close();
}
