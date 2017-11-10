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

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/optional.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/stable_partition.hpp>

#include "cluster/map/map_results.hpp"
#include "cluster/new_cluster_data.hpp"
#include "cluster/old_cluster_data.hpp"
#include "cluster/options/spec/clust_mapping_spec.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/sort_uniq_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/type_aliases.hpp"

using namespace cath;
using namespace cath::clust::detail;
using namespace cath::clust;
using namespace cath::common;
using namespace cath::seq;

using boost::adaptors::filtered;
using boost::algorithm::any_of;
using boost::none;
using boost::range::stable_partition;
using std::less;
using std::max;
using std::string;

/// \brief Map old clusters to new clusters
///
/// \todo Consider trying to break this long function up further
map_results cath::clust::map_clusters(const old_cluster_data_opt &arg_old_clusters,     ///< The old clusters
                                      const new_cluster_data     &arg_new_clusters,     ///< The new clusters
                                      const clust_mapping_spec   &arg_mapping_spec,     ///< The specification for the mapping
                                      const ostream_ref_opt      &arg_domain_out_stream ///< An optional stream to which individual domain mappings should be printed
                                      ) {
	const string dom_map_result_tag = "DOMAIN-MAP-RESULT";

	// Grab the number of new clusters
	const size_t num_new_clusters = get_num_clusters( arg_new_clusters );

	if ( arg_domain_out_stream ) {
		arg_domain_out_stream->get() << "# Columns: tag(" << dom_map_result_tag << ") old_domain cluster_of_old_domain best_overlap_pc (best_new_domain) (cluster_of_best_new_domain) # [where last two columns absent if no new domains on the sequence]\n";
	}

	// Prepare some data structures
	doub_vec           highest_old_clust_overlap_fractions;
	overlap_frac_distn highest_old_dom_overlap_fractions;
	size_t             num_with_nothing_on_parent = 0;
	potential_map_vec  chosen_maps;
	potential_map_vec  potential_maps;
	size_vec           num_mapped_by_new_cluster;
	size_vec           num_mapped_by_old_cluster;

	// If there are old clusters perform a mapping
	if ( arg_old_clusters ) {

		// Grab the number of old clusters and initialise num_mapped_by_new_cluster & num_mapped_by_old_cluster
		const size_t num_old_clusters = arg_old_clusters->size();
		num_mapped_by_new_cluster.resize( num_new_clusters, 0 );
		num_mapped_by_old_cluster.resize( num_old_clusters, 0 );

		// For each old cluster
		for (const size_t &old_cluster_idx : indices( num_old_clusters ) ) {
			const cluster_domains &old_cluster = ( *arg_old_clusters ) [ old_cluster_idx ];

			// Initialise a store of the number of equivalents to this old cluster for each of the new clusters
			size_vec new_clust_equivs( num_new_clusters, 0 );

			// Prepare a function for recording a new domain mapping
			const auto record_mapping_fn = [&] (const size_t &x) {
				++( num_mapped_by_old_cluster[ old_cluster_idx ] );
				++( num_mapped_by_new_cluster[ x               ] );
				++( new_clust_equivs         [ x               ] );
				
			};

			// Loop over the sequences in the old cluster
			//
			// \TODO Come C++17 and structured bindings, use here
			for (const seq_id_and_domain_cluster_ids_pair &old_seq_data : old_cluster) {
				const id_of_string::id_type &seq_id              = old_seq_data.seq_id;
				const domain_cluster_ids    &old_dom_cluster_ids = old_seq_data.dom_cluster_ids;

				// If the new clusters have no entries on the sequence then record overlaps of 0 for all the old domains
				if ( ! has_domain_cluster_ids_of_seq_id( arg_new_clusters, seq_id ) ) {
					num_with_nothing_on_parent += old_dom_cluster_ids.size();
					highest_old_dom_overlap_fractions.add_overlap_fraction( 0.0, old_dom_cluster_ids.size() );

					// If arg_domain_out_stream, print out the name of any old seq ID for which *nothing* could be found
					if ( arg_domain_out_stream ) {
						for (const auto &old_dom_clust_id : old_dom_cluster_ids) {
							arg_domain_out_stream->get()
								<< dom_map_result_tag
								<< " "
								<< arg_old_clusters->get_id_of_seq_name().get_name_of_id( seq_id )
								<< get_segments_suffix_string( old_dom_clust_id.segments )
								<< " "
								<< get_name_of_cluster_of_id( *arg_old_clusters, old_dom_clust_id.cluster_id )
								<< " 0\n";
						}
					}
				}

				// Otherwise, analyse the new clusters' entries on the sequence
				else {
					const domain_cluster_ids &new_dom_clust_ids = get_domain_cluster_ids_of_seq_id( arg_new_clusters, seq_id );

					/// Loop over the old entries on the sequence
					// \TODO Come C++17 and structured bindings, use here
					for (const domain_cluster_id &old_dom_clust_id : old_dom_cluster_ids) {
						const seq_seg_run_opt &old_segments_opt = old_dom_clust_id.segments;

						// If the entry doesn't have segments then...
						if ( ! old_segments_opt ) {

							// Check that old and new both have exactly one entry on this sequence and that the new entry doesn't have segments
							if ( old_dom_cluster_ids.size() != 1 || new_dom_clust_ids.size() != 1 || front( new_dom_clust_ids ).segments ) {
								BOOST_THROW_EXCEPTION(invalid_argument_exception(
									"Inconsistent whole-chain-domain on seq "
									+ arg_old_clusters->get_id_of_seq_name().get_name_of_id( seq_id )
								));
							}

							const auto &equiv_new = front( new_dom_clust_ids );

							// If arg_domain_out_stream, print out the name of any old seq ID for which *nothing* could be found
							if ( arg_domain_out_stream ) {
								arg_domain_out_stream->get()
									<< dom_map_result_tag
									<< " "
									<< arg_old_clusters->get_id_of_seq_name().get_name_of_id( seq_id )
									<< " "
									<< get_name_of_cluster_of_id( *arg_old_clusters, old_dom_clust_id.cluster_id )
									<< " 100 "
									<< arg_old_clusters->get_id_of_seq_name().get_name_of_id( seq_id )
									<< " "
									<< get_name_of_cluster_of_id( arg_new_clusters, equiv_new.cluster_id )
									<< "\n";
							}

							// ...and record that the entries are equivalent 
							record_mapping_fn( equiv_new.cluster_id );
							highest_old_dom_overlap_fractions.add_overlap_fraction( 1.0 );
							continue;
						}

#ifndef NDEBUG
						// Check the old cluster ID is consistent
						if ( old_cluster_idx != old_dom_clust_id.cluster_id ) {
							BOOST_THROW_EXCEPTION(out_of_range_exception("Internal inconsistency detected in old cluster IDs"));
						}
#endif

						// Check that the new_dom_clust_ids isn't empty
						if ( new_dom_clust_ids.empty() ) {
							BOOST_THROW_EXCEPTION(out_of_range_exception("Empty new domain cluster IDs detected"));
						}

						// Check that all of the new_dom_clust_ids have segments
						if ( any_of( new_dom_clust_ids, [&] (const domain_cluster_id &x) { return ! x.segments; } ) ) {
							BOOST_THROW_EXCEPTION(invalid_argument_exception(
								"Inconsistent whole-chain-domain on seq "
								+ arg_old_clusters->get_id_of_seq_name().get_name_of_id( seq_id )
							));
						}

						// Make a closure for calculating the domain overlap with *old_segments_opt
						const auto get_dom_ol_fn = [&] (const domain_cluster_id &x) {
							return fraction_overlap_over_longer( *old_segments_opt, *x.segments );
						};

						// Find the new entry that maps to the old domain best (and its overlap)
						const auto &new_with_best_ol_over_longer = *max_proj_element(
							new_dom_clust_ids,
							less<>{},
							[&] (const domain_cluster_id &x) { return get_dom_ol_fn( x ); }
						);
						const double best_ol = get_dom_ol_fn( new_with_best_ol_over_longer );

						// Record the best overlap
						highest_old_dom_overlap_fractions.add_overlap_fraction( best_ol );

						if ( arg_domain_out_stream ) {
							arg_domain_out_stream->get()
								<< dom_map_result_tag
								<< " "
								<< arg_old_clusters->get_id_of_seq_name().get_name_of_id( seq_id )
								<< get_segments_suffix_string( old_segments_opt )
								<< " "
								<< get_name_of_cluster_of_id( *arg_old_clusters, old_dom_clust_id.cluster_id )
								<< " "
								<< ( 100.0 * best_ol )
								<< " "
								<< arg_old_clusters->get_id_of_seq_name().get_name_of_id( seq_id )
								<< get_segments_suffix_string( new_with_best_ol_over_longer.segments )
								<< " "
								<< get_name_of_cluster_of_id( arg_new_clusters, new_with_best_ol_over_longer.cluster_id )
								<< "\n";
						}

						// If the best entry is above the threshold, record it
						// Note: it is important that this remains a strict >
						// TODO: add a test that demonstrably checks for that
						if ( best_ol > arg_mapping_spec.get_min_equiv_dom_ol() ) {
							record_mapping_fn( new_with_best_ol_over_longer.cluster_id );
						}
					}
				}
			}

			// At the end of the old cluster, store any potential new maps
			for (const size_t &new_cluster_idx : indices( num_new_clusters ) ) {
				if ( new_clust_equivs[ new_cluster_idx ] > 0 ) {
					potential_maps.emplace_back(
						old_cluster_idx,
						new_cluster_idx,
						new_clust_equivs[ new_cluster_idx ]
					);
				}
			}
		}

		// Calculate the highest overlap for each of the old clusters
		highest_old_clust_overlap_fractions.resize( num_old_clusters, 0.0 );
		for (const auto &pot_map : potential_maps) {
			const size_t &old_cluster_idx = pot_map.old_cluster_idx;
			const size_t  num_in_old      = num_entries( ( *arg_old_clusters ) [ old_cluster_idx ] );
			highest_old_clust_overlap_fractions[ old_cluster_idx ] = max(
				highest_old_clust_overlap_fractions[ old_cluster_idx ],
				debug_numeric_cast<double>( pot_map.num_mapped ) / debug_numeric_cast<double>( num_in_old )
			);
		}

		// Stably partition the potential maps to put all the chosen ones at the end and grab the partition point
		const auto partition_point_itr = stable_partition(
			potential_maps,
			[&] (const potential_map &pot_map) {
				const size_t &num_in_old        = num_entries( ( *arg_old_clusters ) [ pot_map.old_cluster_idx ] );
				const size_t &num_in_new        = get_size_of_cluster_of_id( arg_new_clusters, pot_map.new_cluster_idx );
				const size_t &num_mapped        = pot_map.num_mapped;
				const size_t &num_mapped_in_new = num_mapped_by_new_cluster[ pot_map.new_cluster_idx ];

				const double frac_of_old        = debug_numeric_cast<double>( num_mapped ) / debug_numeric_cast<double>( num_in_old        );
				const double frac_of_new        = debug_numeric_cast<double>( num_mapped ) / debug_numeric_cast<double>( num_in_new        );
				const double frac_of_mapped_new = debug_numeric_cast<double>( num_mapped ) / debug_numeric_cast<double>( num_mapped_in_new );

				return ( frac_of_old <= arg_mapping_spec.get_min_equiv_clust_ol()
					||
					frac_of_mapped_new <= clust_mapping_spec::MIN_EQUIV_FRAC_OF_NEW_CLUST_EQUIVS
					||
					frac_of_new        <= clust_mapping_spec::MIN_EQUIV_FRAC_OF_NEW_CLUST
				);
			}
		);

		// Move the chosen maps from the end of potential_maps into chosen_maps
		chosen_maps.reserve( debug_numeric_cast<size_t>( distance( partition_point_itr, end( potential_maps ) ) ) );
		chosen_maps.assign(
			make_move_iterator( partition_point_itr   ),
			make_move_iterator( end( potential_maps ) )
		);
		potential_maps.erase( partition_point_itr, common::cend( potential_maps ) );
	}

	// Return the results of the mapping
	return {
		chosen_maps,
		potential_maps,
		get_info_ordered_indices_of_unmapped_new_clusters(
			chosen_maps,
			arg_new_clusters
		),
		num_mapped_by_new_cluster,
		num_mapped_by_old_cluster,
		num_with_nothing_on_parent,
		std::move( highest_old_dom_overlap_fractions ),
		build_overlap_frac_distn_from_overlap_fractions( highest_old_clust_overlap_fractions ),
		arg_mapping_spec
	};
}

/// \brief Get the indices of the specified new clusters not mentioned in the specified chosen maps,
///        ordered by the info of the corresponding clusters
size_vec cath::clust::get_info_ordered_indices_of_unmapped_new_clusters(const potential_map_vec &arg_chosen_maps, ///< The chosen maps
                                                                        const new_cluster_data  &arg_new_clusters ///< The new clusters
                                                                        ) {
	const size_t num_new_clusters = get_num_clusters( arg_new_clusters );

	// Build a deque indicating which of the new clusters are mapped
	bool_deq mapped_new_clusters( num_new_clusters, false );
	for (const potential_map &x : arg_chosen_maps) {
		mapped_new_clusters[ x.new_cluster_idx ] = true;
	}

	// Build a size_vec of the indices of the unmapped new clusters, sorted using cluster_info
	return sort_build<size_vec>(
		indices( num_new_clusters )
			| filtered( [&] (const size_t &x) { return ! mapped_new_clusters[ x ]; } ),
		[&] (const size_t &x, const size_t &y) {
			return (
				get_info_of_cluster_of_id( arg_new_clusters, x )
				<
				get_info_of_cluster_of_id( arg_new_clusters, y )
			);
		}
	);
}
