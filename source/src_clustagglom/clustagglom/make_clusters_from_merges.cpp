/// \file
/// \brief The make_clusters_from_merges class definitions

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

#include "make_clusters_from_merges.hpp"

#include <boost/optional.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/range/adaptor/filtered.hpp>

#include "clustagglom/detail/node_info.hpp"
#include "clustagglom/hierarchy.hpp"
#include "clustagglom/merge.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/debug_numeric_cast.hpp"

#include <iostream>

using namespace cath::clust;
using namespace cath::clust::detail;
using namespace cath::common;

using boost::adaptors::filtered;
using boost::adaptors::reversed;
using boost::none;
using boost::remove_erase_if;
using boost::sub_range;
using std::cref;
using std::max;
using std::move;

/// \brief Make a vector of iterators into the specified vector of merges that demarcate the
///        boundaries of the each of the regions up to (and equal to) each of the successively
///        looser cutoffs.
///
/// This starts at begin( arg_merges ) and then has one iterator for each cutoff after that.
/// (No end is require because merges after the last cutoff aren't considered)
///
/// Merges that have a dissimilarity exactly equalling a cutoff will be included in the region before
/// not excluded and pushed out or into the region after (not the use of `upper_bound()` not `lower_bound()`)
///
/// '''upper_bound'''
///
/// \pre `is_sorted( arg_cutoffs );`
merge_vec_citr_vec cath::clust::detail::calc_merge_cutoff_boundaries(const merge_vec    &arg_merges, ///< The ordered list of merges to determine the clusters
                                                                     const strength_vec &arg_cutoffs ///< The cutoffs at which the clusters should be drawn, in ascending order
                                                                     ) {
	merge_vec_citr_vec merge_cutoff_boundaries;
	merge_cutoff_boundaries.reserve  ( arg_cutoffs.size() + 1       );
	merge_cutoff_boundaries.push_back( common::cbegin( arg_merges ) );
	for (const strength &cutoff : arg_cutoffs) {
		merge_cutoff_boundaries.push_back( upper_bound(
			merge_cutoff_boundaries.back(),
			common::cend( arg_merges ),
			cutoff,
			[] (const strength &x, const merge &y) { return x < y.dissim; }
		) );
	}

	return merge_cutoff_boundaries;
}

/// \brief Make clusters (represented in a hierarchy object) at the specified cutoffs from
///        the specified merges over the specified number of items
///
/// The merge_vec must follow the following rules:
///  * the indices of leaf nodes are all less than arg_num_entities
///  * the indices of all non-leaf nodes are all >= arg_num_entities and < 2*arg_num_entities
///  * the index of a parent node is > than both indices of the children nodes
///
/// \pre `is_sorted( arg_cutoffs )`
hierarchy cath::clust::make_clusters_from_merges(const merge_vec    &arg_merges,       ///< The ordered list of merges to determine the clusters
                                                 const item_idx     &arg_num_entities, ///< The number of items being clustered
                                                 const strength_vec &arg_cutoffs       ///< The cutoffs at which the clusters should be drawn, in ascending order
                                                 ) {
	using std::to_string;

	// The layers that will be populated
	hierarchy_layer_vec hier_layers;
	const size_t &num_merges   = arg_merges.size();
	const size_t  num_entities = arg_num_entities;
	item_opt_vec root_parent_of_node( num_entities + max( num_entities, num_merges ) );

	{
		const auto merge_cutoff_boundaries = calc_merge_cutoff_boundaries( arg_merges, arg_cutoffs );

		// Prepare root_parent_of_node for calculation
		node_info_opt_vec info_of_cutoff_root( num_entities + max( num_entities, num_merges ) );
		item_vec          new_clust_indices;
		item_vec          unused_prev_clust_indices;

		// Loop over each over the regions before each of the cutoffs
		for (const size_t &cutoff_layer_ctr : indices( arg_cutoffs.size() ) ) {

			// Prepare a sub_range for this region of the merges
			const sub_range<const merge_vec> cutoff_subrange_of_merges{
				merge_cutoff_boundaries[ cutoff_layer_ctr     ],
				merge_cutoff_boundaries[ cutoff_layer_ctr + 1 ]
			};

			// Looping over merges (normal order, ie best-to-worst)...
			// make the slots in root_parent_of_node for each of the child nodes contain the index of their parent
			for (const auto &x : cutoff_subrange_of_merges) {
				for (const auto &child_node : { cref( x.node_a ), cref( x.node_b ) } ) {
					root_parent_of_node[ child_node ] = x.merge_node;
				}
			}

			// Looping over merges, (*reverse order, ie worst-to-best*)...
			for (const auto &x : cutoff_subrange_of_merges | reversed ) {

				// If this is a root merge within this cutoff regions
				// (because it doesn't have an entry in root_parent_of_node for its new node):
				//    * Create and populate the node_info slot in info_of_cutoff_root for the root node
				//    * Add the index of the root node to new_clust_indices
				//    * Add a new group to layer
				if ( ! static_cast<bool>( root_parent_of_node[ x.merge_node ] ) ) {
					info_of_cutoff_root[ x.merge_node ] = node_info( new_clust_indices.size(), cutoff_layer_ctr );
					new_clust_indices.push_back( x.merge_node );
				}

				// Else propagate the value in the slot in root_parent_of_node for the parent node to the slots
				// for the child nodes
				//
				// Because of the ordering, this puts the index of each root node (within this cutoff region)
				// in all of the slots in root_parent_of_node corresponding to their children
				else {
					root_parent_of_node[ x.node_a ] = *root_parent_of_node[ x.merge_node ];
					root_parent_of_node[ x.node_b ] = *root_parent_of_node[ x.merge_node ];
				}
			}

			// Looping over merges, (normal order, ie best-to-worst)...
			hierarchy_layer layer( new_clust_indices.size() );
			for (const auto &x : cutoff_subrange_of_merges) {

				// Get the index of the hierarcy_group for this merge's root node
				const size_t &group_index = info_of_cutoff_root[ *root_parent_of_node[ x.node_a ] ]->locn;

				// For each of the two child nodes of the merge
				for (const auto &child_node : { cref( x.node_a ), cref( x.node_b ) } ) {
					const node_info_opt &child_info_opt = info_of_cutoff_root[ child_node ];

					// If this child-node is a leaf node add an entry for it in the group
					if ( child_node < arg_num_entities ) {
						layer[ group_index ].emplace_back( hierarchy_ref::ENTRY,   child_node );
					}
					// Else if it's a root node of the previous cutoff region...
					// *** WHY IS THIS *THE* PREVIOUS NODE? ***
					else if ( static_cast<bool>( child_info_opt ) ) {
						if ( child_info_opt->layer + 1 == cutoff_layer_ctr ) {
							layer[ group_index ].emplace_back( hierarchy_ref::CLUSTER, child_info_opt->locn );
						}
					}
				}
			}

			// Looping over merges, (normal order, ie best-to-worst)...
			// wipe the info for any previous root nodes merged in this cutoff region
			for (const auto &x : cutoff_subrange_of_merges) {
				for (const auto &child_node : { cref( x.node_a ), cref( x.node_b ) } ) {
					if ( static_cast<bool>( info_of_cutoff_root[ child_node ] ) ) {
						info_of_cutoff_root[ child_node ] = none;
					}
				}
			}

			// Remove any of the unused_prev_clust_indices that have been used in this cutoff region
			remove_erase_if(
				unused_prev_clust_indices,
				[&] (const size_t &x) {
					return ! static_cast<bool>( info_of_cutoff_root[ x ] );
				}
			);

			// Propagate any unused_prev_clusters through a new separate group
			for (const item_idx &previous_clust_num : unused_prev_clust_indices) {
				node_info &unused_prev_info = *info_of_cutoff_root[ previous_clust_num ];
				const size_t prev_locn = unused_prev_info.locn;

				unused_prev_info.layer = cutoff_layer_ctr;
				unused_prev_info.locn  = layer.size();

				layer.emplace_back( move( hierarchy_group{}.emplace_back( hierarchy_ref::CLUSTER, prev_locn ) ) );
			}

			// Add the numbers of this cutoff-region's new clusters
			unused_prev_clust_indices.insert(
				common::cend  ( unused_prev_clust_indices ),
				common::cbegin( new_clust_indices         ),
				common::cend  ( new_clust_indices         )
			);

			// Clear the new_clust_indices for reuse in the next cutoff-region
			new_clust_indices.clear();

			// Add this to the list of layers
			hier_layers.push_back( move( layer ) );
		}
	}

	// Return the result of adding a root layer, reversing and using to construct a hierarchy
	return make_hierarchy_from_reversed_without_root_layer(
		hier_layers,
		indices( static_cast<size_t>( arg_num_entities ) )
			| filtered( [&] (const size_t &x) { return ! root_parent_of_node[ x ]; } )
	);
}

/// \brief Make clusters (represented in a hierarchy object) at the specified cutoffs from
///        the specified merges and then sort it based on the specified sorting indices for the items
///
/// \pre `is_sorted( arg_cutoffs )`
hierarchy cath::clust::make_clusters_from_merges_and_sort(const merge_vec    &arg_merges,          ///< The ordered list of merges to determine the clusters
                                                          const size_vec     &arg_sorting_indices, ///< Values corresponding to each of the entries such that one value is less than another if the corresponding entry should be sorted before the other
                                                          const strength_vec &arg_cutoffs          ///< The cutoffs at which the clusters should be drawn, in ascending order
                                                          ) {
	return sort_hierarchy_copy(
		make_clusters_from_merges(
			arg_merges,
			debug_numeric_cast<item_idx>( arg_sorting_indices.size() ),
			arg_cutoffs
		),
		arg_sorting_indices
	);

}
