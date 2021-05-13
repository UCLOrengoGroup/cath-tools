/// \file
/// \brief The calc_complete_linkage_merge_list class definitions

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

#include "calc_complete_linkage_merge_list.hpp"

#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/range/algorithm/partition.hpp>
#include <boost/range/algorithm/upper_bound.hpp>

#include "cath/clustagglom/detail/clust_id_pot.hpp"
#include "cath/clustagglom/links.hpp"
#include "cath/clustagglom/merge.hpp"
#include "cath/common/boost_addenda/range/stable_sort_proj.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/optional/make_optional_if.hpp"

using namespace ::cath::clust;
using namespace ::cath::common;

using ::boost::range::for_each;
using ::boost::range::min_element;
using ::boost::range::partition;
using ::boost::range::upper_bound;
using ::std::make_optional;
using ::std::max;
using ::std::min;
using ::std::tie;

/// \brief Calculate the ordered sequence of merges to be conducted by complete-linkage clustering
///        given the specified links and item ordering
///
/// Note that the result is just a sequence of merges (in descending order of quality),
/// not (yet) a list of clusters. Use make_clusters_from_merges() on this output to get clusters.
///
/// \TODO Consider adding a parameter of the worst tolerable dissimilarity
///       so that effort needn't be wasted calculating merges past that point
///
/// \relates links
merge_vec cath::clust::calc_complete_linkage_merge_list(links             prm_links,        ///< The links to analyse
                                                        const size_vec   &prm_sort_indices, ///< The ranks of the items (ie a 0 should appear in the index corresponding to that of the most preferred item)
                                                        const strength   &prm_max_dissim    ///< The maximum dissimilarity at which merges may still happen
                                                        ) {
	using ::std::to_string;

	const size_t num_entities = prm_sort_indices.size();
	// std::cerr << "num_entities is : " << num_entities << "\n";

	// const auto cluster_start_time = high_resolution_clock::now();

	merge_vec results;

	detail::clust_id_pot clust_ids( num_entities );
	size_t num_clusts = num_entities;

	if ( ! prm_links.empty() ) {

		size_vec sorted_indices{ prm_sort_indices };
		size_vec sizes( prm_links.size(), 1_z );
		item_vec chain;

		while ( num_clusts > 1 ) {
			const bool start_new_chain = ( chain.size() < 4_z );
			item_idx a, b;
			if ( start_new_chain ) {
				a = clust_ids.get_jumbled_nth_index( 0 );
				b = clust_ids.get_jumbled_nth_index( 1 );
				chain.assign( 1, a );
			}
			else {
				a = chain[ chain.size() - 4_z ];
				b = chain[ chain.size() - 3_z ];

				// The round-trip static_cast to uint32_t and then back to size_t
				// for the call to resize() shouldn't be necessary but otherwise GCC 7.2
				// gives a false-positive warning:
				//
				//     In function ‘cath::clust::merge_vec: cath::clust::calc_complete_linkage_merge_list(cath::clust::links, const size_vec&, const strength&)’:
				//     cc1plus: error: ‘void* __builtin_memset(void*, int, long unsigned int)’: specified size 18446744073709551604 exceeds maximum object size 9223372036854775807 [-Werror=stringop-overflow=]
				//
				// The warning isn't very descriptive and doesn't give a location but
				// appears to be about this resize() argument being too big. Presumably, it's deducing
				// that the value might wrap below 0 to a huge value (because it's an unsigned integer type)
				// but we know that, within this else clause, chain.size() is big enough to prevent that.
				//
				// I've reported this issue here:
				//   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=83239
				chain.resize( static_cast<uint32_t>( chain.size() - 3_z ) );
			}

			strength dist;
			do {
				auto &a_links = prm_links[ a ];

				a_links.erase(
					partition(
						a_links,
						[&] (const link &x) { return clust_ids.has_index( x.node ); }
					),
					common::cend( a_links )
				);
				const auto min_itr = min_element(
					a_links,
					[&] (const link &x, const link &y) {
						const char x_is_b_score = ( ( x.node == b ) ? 0 : 1 );
						const char y_is_b_score = ( ( y.node == b ) ? 0 : 1 );
						return (
							tie( x.dissim, sorted_indices[ x.node ], x_is_b_score )
							<
							tie( y.dissim, sorted_indices[ y.node ], y_is_b_score )
						);
					}
				);

				b = a;
				if ( min_itr == common::cend( a_links ) ) {
					// BOOST_THROW_EXCEPTION(out_of_range_exception("argh"));
					a    = clust_ids.get_min_value_excluding_spec( a );
					dist = std::numeric_limits< decltype( dist ) >::infinity();
				}
				else {
					a    = min_itr->node;
					dist = min_itr->dissim;
				}
				chain.push_back( a );
			} while ( chain.size() < 3 || a != chain[ chain.size() - 3 ] );

			clust_ids.remove_index( a )
			         .remove_index( b );

			const item_idx &new_label = clust_ids.add_new_index();
			if ( sizes.size() < new_label + 1 ) {
				sizes.resize( new_label + 1 );
			}

			results.emplace_back( a, b, new_label, dist );
			--num_clusts;


			sizes[ new_label ] = sizes[ a ] + sizes[ b ];

			// Update sorted_indices
			if ( new_label != sorted_indices.size() ) {
				BOOST_THROW_EXCEPTION(out_of_range_exception("ARGH")); /// \TODO Remove this if it doesn't fire for a while
			}
			sorted_indices.push_back( min( sorted_indices[ a ], sorted_indices[ b ] ) );

			prm_links.merge(
				a,
				b,
				new_label,
				// !!!!! At the moment, this function is only called where both links exist - OK for complete-linkage but not others
				[&] (const item_idx &x, ///< The index of the target cluster
				     const strength &y, ///< The dissimilarity that the first  cluster had to the target cluster
				     const strength &z  ///< The dissimilarity that the second cluster had to the target cluster
				     ) {
					return if_then_optional(
						( clust_ids.has_index( x ) ),
						max( y, z )
					);
				}
			);
		}

		// Stable-sort the results
		stable_sort_proj(
			results,
			std::less<>{},
			[] (const merge &x) { return x.dissim; }
		);

		// Remove any merges that went beyond prm_max_dissim
		results.erase(
			upper_bound(
				results,
				prm_max_dissim,
				[] (const strength &max_dissim, const merge &x) { return max_dissim < x.dissim; }
			),
			common::cend( results )
		);
	}

	// Ensure that each merge has the lower node ID first, swapping as necessary
	//
	// (Should take very little effort and helps to make merge list more reproducible)
	for_each(
		results,
		[] (merge &x) {
			if ( x.node_a > x.node_b ) {
				std::swap( x.node_a, x.node_b );
			}
		}
	);

	// Return the results
	return results;
}


/// \brief Calculate the ordered sequence of merges to be conducted by complete-linkage clustering
///        given the specified links
///
/// Note that the result is just a sequence of merges (in descending order of quality),
/// not (yet) a list of clusters. Use make_clusters_from_merges() on this output to get clusters.
///
/// \TODO Consider adding a parameter of the worst tolerable dissimilarity
///       so that effort needn't be wasted calculating merges past that point
///
/// Since no preferred ranking of the items is specified, any ambiguities will
/// be resolved by preferring the items in descending order
merge_vec cath::clust::calc_complete_linkage_merge_list(links             prm_links,     ///< The links to analyse
                                                        const size_t     &prm_size,      ///< The number of items to be merged
                                                        const strength   &prm_max_dissim ///< The maximum dissimilarity at which merges may still happen
                                                        ) {
	return calc_complete_linkage_merge_list(
		std::move( prm_links ),
		copy_build<size_vec>( indices( prm_size ) ),
		prm_max_dissim
	);
}

/// \brief Calculate the ordered sequence of merges to be conducted by complete-linkage clustering
///        given the specified links and item ordering
///
/// Note that the result is just a sequence of merges (in descending order of quality),
/// not (yet) a list of clusters. Use make_clusters_from_merges() on this output to get clusters.
///
/// \TODO Consider adding a parameter of the worst tolerable dissimilarity
///       so that effort needn't be wasted calculating merges past that point
merge_vec cath::clust::calc_complete_linkage_merge_list(const item_item_strength_tpl_vec &prm_links,        ///< The links to analyse
                                                        const size_vec                   &prm_sort_indices, ///< The ranks of the items (ie a 0 should appear in the index corresponding to that of the most preferred item)
                                                        const strength                   &prm_max_dissim    ///< The maximum dissimilarity at which merges may still happen
                                                        ) {
	return calc_complete_linkage_merge_list(
		make_links( prm_links ),
		prm_sort_indices,
		prm_max_dissim
	);
}

/// \brief Calculate the ordered sequence of merges to be conducted by complete-linkage clustering
///        given the specified links
///
/// Note that the result is just a sequence of merges (in descending order of quality),
/// not (yet) a list of clusters. Use make_clusters_from_merges() on this output to get clusters.
///
/// \TODO Consider adding a parameter of the worst tolerable dissimilarity
///       so that effort needn't be wasted calculating merges past that point
///
/// Since no preferred ranking of the items is specified, any ambiguities will
/// be resolved by preferring the items in descending order
merge_vec cath::clust::calc_complete_linkage_merge_list(const item_item_strength_tpl_vec &prm_links,     ///< The links to analyse
                                                        const size_t                     &prm_size,      ///< The number of items to be merged
                                                        const strength                   &prm_max_dissim ///< The maximum dissimilarity at which merges may still happen
                                                        ) {
	return calc_complete_linkage_merge_list(
		make_links( prm_links ),
		prm_size,
		prm_max_dissim
	);
}
