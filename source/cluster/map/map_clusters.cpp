/// \map
/// \brief The map_clusters class definitions

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

#include "map_clusters.hpp"

#include <boost/optional.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/stable_partition.hpp>
#include <boost/range/irange.hpp>

#include "cluster/map/map_results.hpp"
#include "cluster/new_cluster_data.hpp"
#include "cluster/old_cluster_data.hpp"
#include "cluster/options/spec/clust_mapping_spec.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/sort_uniq_build.hpp"
#include "common/type_aliases.hpp"

using namespace cath::clust::detail;
using namespace cath::clust;
using namespace cath::common;
using namespace cath::seq;

using boost::adaptors::filtered;
using boost::find_if;
using boost::irange;
using boost::none;
using boost::range::stable_partition;
using boost::sub_range;
using std::deque;
using std::vector;

/// \brief Map old clusters to new clusters
map_results cath::clust::map_clusters(const old_cluster_data_opt &arg_old_clusters, ///< The old clusters
                                      const new_cluster_data     &arg_new_clusters, ///< The new clusters
                                      const clust_mapping_spec   &arg_mapping_spec  ///< The specification for the mapping
                                      ) {
	const size_t num_new_clusters = get_num_clusters( arg_new_clusters );

	const double &min_equiv_dom_ol   = arg_mapping_spec.get_min_equiv_dom_ol();
	const double &min_equiv_clust_ol = arg_mapping_spec.get_min_equiv_clust_ol();

	// size_opt_vec old_clust_id_of_new_clust_id;
	size_vec num_mapped_by_new_cluster;
	size_vec num_mapped_by_old_cluster;
	// old_clust_id_of_new_clust_id.resize( num_new_clusters, none );
	num_mapped_by_new_cluster.resize( num_new_clusters, 0 );

	potential_map_vec chosen_maps;
	potential_map_vec potential_maps;

	doub_vec domain_mapping_fractions;
	doub_vec cluster_mapping_fractions;

	if ( arg_old_clusters ) {
		const size_t num_old_clusters = arg_old_clusters->size();
		num_mapped_by_old_cluster.resize( num_old_clusters, 0 );

		for (const size_t &old_cluster_idx : irange( 0_z, num_old_clusters ) ) {
			const cluster_domains &old_cluster = ( *arg_old_clusters ) [ old_cluster_idx ];
			// const auto dom_require;
			// min_equiv_dom_ol;
			// min_equiv_dom_ol;

			size_vec new_clust_equivs( num_new_clusters, 0 );

			const auto record_mapping_fn = [&] (const size_t &x) {
				++( num_mapped_by_old_cluster[ old_cluster_idx ] );
				++( num_mapped_by_new_cluster[ x               ] );
				++( new_clust_equivs         [ x               ] );
			};

			// \TODO Come C++17 and structured bindings, use here
			for (const seq_id_and_domain_cluster_ids_pair &old_seq_data : old_cluster) {
				const id_of_string::id_type &seq_id              = old_seq_data.seq_id;
				const domain_cluster_ids    &old_dom_cluster_ids = old_seq_data.dom_cluster_ids;

				if ( has_domain_cluster_ids_of_seq_id( arg_new_clusters, seq_id ) ) {
					const domain_cluster_ids &new_dom_clust_ids = get_domain_cluster_ids_of_seq_id( arg_new_clusters, seq_id );

					// \TODO Come C++17 and structured bindings, use here
					for (const domain_cluster_id &old_dom_clust_id : old_dom_cluster_ids) {
						const seq_seg_run_opt &old_segments_opt = old_dom_clust_id.segments;


						if ( ! old_segments_opt ) {
							if ( old_dom_cluster_ids.size() != 1 || new_dom_clust_ids.size() != 1 || front( new_dom_clust_ids ).segments ) {
								BOOST_THROW_EXCEPTION(invalid_argument_exception("Inconsistent whole-chain-domain"));
							}
							record_mapping_fn( front( new_dom_clust_ids ).cluster_id );
							continue;
						}

#ifndef NDEBUG
						if ( old_cluster_idx != old_dom_clust_id.cluster_id ) {
							BOOST_THROW_EXCEPTION(out_of_range_exception("Internal inconsistency detected in old cluster IDs"));
						}
#endif
						const auto find_new_itr = find_if(
							new_dom_clust_ids,
							[&] (const domain_cluster_id &x) {
								if ( ! x.segments ) {
									BOOST_THROW_EXCEPTION(invalid_argument_exception("Inconsistent whole-chain-domain"));
								}
								return (
									fraction_overlap_over_longer( *old_segments_opt, *x.segments )
									>=
									min_equiv_dom_ol
								);
							}
						);
						if ( find_new_itr != common::cend( new_dom_clust_ids ) ) {
							record_mapping_fn( find_new_itr->cluster_id );
						}
						// fraction_overlap_over_longer;
						// new_dom_clust_ids;
						// old_segments_opt;
					}

				}
			}

			for (const size_t &new_cluster_idx : irange( 0_z, num_new_clusters ) ) {
				if ( new_clust_equivs[ new_cluster_idx ] > 0 ) {
					chosen_maps.emplace_back(
						old_cluster_idx,
						new_cluster_idx,
						new_clust_equivs[ new_cluster_idx ]
					);
				}
			}

		}

		const auto partition_point_itr = stable_partition(
			chosen_maps,
			[&] (const potential_map &pot_map) {
				const size_t &old_cluster_idx   = pot_map.old_cluster_idx;
				const size_t &new_cluster_idx   = pot_map.new_cluster_idx;

				// const size_t &num_in_new        = get_size_of_cluster_of_id( arg_new_clusters, new_cluster_idx );
				const size_t &num_in_old        = ( *arg_old_clusters ) [ old_cluster_idx ].size();
				const size_t &num_mapped        = pot_map.num_mapped;
				const size_t &num_mapped_in_new = num_mapped_by_new_cluster[ new_cluster_idx ];
				// const size_t &num_mapped_in_old = num_mapped_by_old_cluster[ old_cluster_idx ];

				const double frac_of_old        = debug_numeric_cast<double>( num_mapped ) / debug_numeric_cast<double>( num_in_old        );
				const double frac_of_mapped_new = debug_numeric_cast<double>( num_mapped ) / debug_numeric_cast<double>( num_mapped_in_new );

				return ( frac_of_old > min_equiv_clust_ol && frac_of_mapped_new > 0.5 );
			}
		);

		sub_range<potential_map_vec> rejected_range{ partition_point_itr, end( chosen_maps ) };
		potential_maps.reserve( rejected_range.size() );
		for (auto &the_map : rejected_range) {
			potential_maps.emplace_back( std::move( the_map ) );
		}
		chosen_maps.erase( partition_point_itr, common::cend( chosen_maps ) );

	}

	deque<bool> mapped_new_clusters( num_new_clusters, false );
	for (const potential_map &x : chosen_maps) {
		mapped_new_clusters[ x.new_cluster_idx ] = true;
	}


	const auto unmapped_new_cluster_indices = sort_build<size_vec>(
		irange( 0_z, num_new_clusters )
			| filtered( [&] (const size_t &x) { return ! mapped_new_clusters[ x ]; } ),
		[&] (const size_t &x, const size_t &y) {
			return (
				get_info_of_cluster_of_id( arg_new_clusters, x )
				<
				get_info_of_cluster_of_id( arg_new_clusters, y )
			);
		}
	);

	return {
		chosen_maps,
		potential_maps,
		unmapped_new_cluster_indices,
		num_mapped_by_new_cluster,
		num_mapped_by_old_cluster,
		domain_mapping_fractions,
		cluster_mapping_fractions,
		arg_mapping_spec
	};
}
