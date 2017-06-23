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

// #include <boost/algorithm/string/join.hpp> // ***** TEMPORARY *****
#include <boost/optional.hpp>
// #include <boost/range/adaptor/transformed.hpp> // ***** TEMPORARY *****
#include <boost/range/irange.hpp>

#include "cluster/new_cluster_data.hpp"
#include "cluster/old_cluster_data.hpp"
#include "cluster/options/spec/clust_mapping_spec.hpp"
#include "common/type_aliases.hpp"

using namespace cath::clust::detail;
using namespace cath::clust;
using namespace cath::common;
using namespace cath::seq;

// using boost::adaptors::transformed; // ***** TEMPORARY *****
// using boost::algorithm::join; // ***** TEMPORARY *****
using boost::find_if;
using boost::irange;
using boost::none;

/// \brief Map old clusters to new clusters
int cath::clust::map_clusters(const old_cluster_data_opt &arg_old_clusters, ///< The old clusters
                              const new_cluster_data     &arg_new_clusters, ///< The new clusters
                              const clust_mapping_spec   &arg_mapping_spec  ///< The specification for the mapping
                              ) {
	const size_t num_new_clusters = get_num_clusters( arg_new_clusters );
	size_opt_vec old_clust_id_of_new_clust_id;
	old_clust_id_of_new_clust_id.resize( num_new_clusters, none );

	const double &min_equiv_dom_ol   = arg_mapping_spec.get_min_equiv_dom_ol();
	// const double &min_equiv_clust_ol = arg_mapping_spec.get_min_equiv_clust_ol();

	if ( arg_old_clusters ) {
		for (const size_t &old_cluster_idx : irange( 0_z, arg_old_clusters->size() ) ) {
			const cluster_domains &old_cluster = ( *arg_old_clusters ) [ old_cluster_idx ];
			// const auto dom_require;
			// min_equiv_dom_ol;
			// min_equiv_dom_ol;

			size_vec new_clust_equivs( num_new_clusters, 0 );

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
							++( new_clust_equivs[ front( new_dom_clust_ids ).cluster_id ] );
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
							++( new_clust_equivs[ find_new_itr->cluster_id ] );
						}
						// fraction_overlap_over_longer;
						// new_dom_clust_ids;
						// old_segments_opt;
					}

				}
			}

			// std::cerr << "For old cluster ";
			// std::cerr << old_cluster_idx;
			// std::cerr << " : ";
			// std::cerr << join(
			// 	new_clust_equivs
			// 		| transformed( [] (const size_t &x) { using std::to_string; return to_string( x ); } ),
			// 	", "
			// );
			// std::cerr << "\n";

			// /// \brief The domain_cluster_ids of the domains on the sequence
			// domain_cluster_ids dom_cluster_ids;

			// For each old cluster member:
			//   Search in the unordered_map to find a domain equiv
			//   If a domain equivalent is found:
			//     Find equivalent domain's cluster
			//     Increment that cluster's counter of domain equivalents
			//     If the counter has reached the target, add the cluster mapping to the list and move to the next cluster

			// arg_mapping_spec;
			// const size_t &target_equivalences = old_cluster;
			// *arg_old_clusters;
		}
	}

	// sort_any_any_unmapped_new_clusters();
	// give_unmapped_new_clusters_a_name()

	return 0;
}
