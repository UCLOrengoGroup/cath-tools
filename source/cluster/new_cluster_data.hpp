/// \file
/// \brief The new_cluster_data class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_NEW_CLUSTER_DATA_H
#define _CATH_TOOLS_SOURCE_CLUSTER_NEW_CLUSTER_DATA_H

#include <boost/range/adaptor/transformed.hpp>

#include "cluster/clusters_info.hpp"
#include "cluster/domain_cluster_ids_by_seq.hpp"
#include "common/boost_addenda/range/accumulate_proj.hpp"

namespace cath {
	namespace clust {

		/// \brief Store the data associated with new "to" clusters
		///
		/// This uses a reference to an external id_of_string for mapping
		/// from sequences names to IDs. This allows the same IDs to be shared
		/// with other data structures.
		class new_cluster_data final {
		private:
			/// \brief Cluster IDs of each of the entries, indexed by their parent entry (sequence)
			domain_cluster_ids_by_seq dom_clust_ids;

			/// \brief Basic information on the cluster
			clusters_info clust_info;

		public:
			explicit new_cluster_data(common::id_of_string &) noexcept;

			/// Prevent construction from an id_of_string rvalue
			new_cluster_data(const common::id_of_string &&) = delete;

			clust_entry_problem add_entry(const boost::string_ref &,
			                              const boost::string_ref &,
			                              const boost::string_ref &,
			                              seq::seq_seg_run_opt);

			const domain_cluster_ids_by_seq & get_dom_clust_ids() const;
			const clusters_info & get_clust_info() const;
		};

		/// \brief Ctor from an id_of_string reference
		inline new_cluster_data::new_cluster_data(common::id_of_string &arg_id_of_str ///< The id_of_string to use to map from sequences names to IDs
		                                          ) noexcept : dom_clust_ids{ arg_id_of_str } {
		}

		/// \brief Add an entry with the specified sequence name and (optional) segments to the cluster with
		///        the specified name
		inline clust_entry_problem new_cluster_data::add_entry(const boost::string_ref &arg_clust_name, ///< The name of the cluster of the entry
		                                                       const boost::string_ref &arg_seq_id,     ///< The name of the sequence within which this entry appears
		                                                       const boost::string_ref &arg_domain_id,  ///< The name of the entry
		                                                       seq::seq_seg_run_opt     arg_segments    ///< The (optional) segments of the entry within the sequence
		                                                       ) {
			try {
				const auto cluster_id = update_info_and_get_id_for_cluster_of_name(
					clust_info,
					arg_clust_name,
					arg_domain_id,
					arg_segments
				);
				return dom_clust_ids.add(
					arg_seq_id,
					domain_cluster_id{ std::move( arg_segments ), cluster_id }
				);
			}
			catch (...) {
				return clust_entry_problem::CLASH;
			}
		}

		/// \brief Get the dom_clust_ids
		inline const domain_cluster_ids_by_seq & new_cluster_data::get_dom_clust_ids() const {
			return dom_clust_ids;
		}

		/// \brief Get the cluster info
		inline const clusters_info & new_cluster_data::get_clust_info() const {
			return clust_info;
		}

		/// \brief Get the number of clusters from the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline size_t get_num_clusters(const new_cluster_data &arg_new_cluster_data ///< The new_cluster_data to query
		                               ) {
			return arg_new_cluster_data.get_clust_info().get_num_clusters();
		}

		/// \brief Get the info of the cluster with the specified ID from the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline const cluster_info & get_info_of_cluster_of_id(const new_cluster_data &arg_new_cluster_data, ///< The new_cluster_data to query
		                                                      const cluster_id_t     &arg_cluster_id        ///< The ID of the cluster of interest
		                                                      ) {
			return arg_new_cluster_data.get_clust_info().get_info_of_cluster_of_id( arg_cluster_id );
		}

		/// \brief Get the size of the cluster with the specified ID from the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline size_t get_size_of_cluster_of_id(const new_cluster_data &arg_new_cluster_data, ///< The new_cluster_data to query
		                                        const cluster_id_t     &arg_cluster_id        ///< The ID of the cluster of interest
		                                        ) {
			return get_info_of_cluster_of_id( arg_new_cluster_data, arg_cluster_id ).get_size();
		}

		/// \brief Get the total number of entries in the clusters of the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline size_t get_num_entries(const new_cluster_data &arg_new_cluster_data ///< The new_cluster_data to query
		                              ) {
			return common::accumulate_proj(
				boost::irange( 0_z, get_num_clusters( arg_new_cluster_data ) ),
				0_z,
				std::plus<>{},
				[&] (const size_t &x) { return get_size_of_cluster_of_id( arg_new_cluster_data, x ); }
			);
		}

		/// \brief Get the name of the cluster with the specified ID in the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline const std::string & get_name_of_cluster_of_id(const new_cluster_data &arg_new_cluster_data, ///< The new_cluster_data to query
		                                                     const cluster_id_t     &arg_cluster_id        ///< The ID of the cluster of interest
		                                                     ) {
			return arg_new_cluster_data.get_clust_info().get_name_of_cluster_of_id( arg_cluster_id );
		}

		/// \brief Get the map from sequence name to sequence ID from the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline const common::id_of_string & get_id_of_seq_name(const new_cluster_data &arg_new_cluster_data ///< The new_cluster_data to query
		                                                       ) {
			return arg_new_cluster_data.get_dom_clust_ids().get_id_of_seq_name();
		}

		/// \brief Get whether there is a non-empty domain_cluster_ids associated with the specified sequence id in the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline bool has_domain_cluster_ids_of_seq_id(const new_cluster_data              &arg_new_cluster_data, ///< The new_cluster_data to query
		                                             const common::id_of_string::id_type &arg_seq_id            ///< The id of the sequence of interest
		                                             ) {
			return has_domain_cluster_ids_of_seq_id( arg_new_cluster_data.get_dom_clust_ids(), arg_seq_id );
		}

		/// \brief Get the domain_cluster_ids associated with the specified sequence id in the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline const domain_cluster_ids & get_domain_cluster_ids_of_seq_id(const new_cluster_data              &arg_new_cluster_data, ///< The new_cluster_data to query
		                                                                   const common::id_of_string::id_type &arg_seq_id            ///< The id of the sequence of interest
		                                                                   ) {
			return arg_new_cluster_data.get_dom_clust_ids()[ arg_seq_id ];
		}

		/// \brief Get whether there is a non-empty domain_cluster_ids associated with the specified sequence name in the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline bool has_domain_cluster_ids_of_seq_name(const new_cluster_data &arg_new_cluster_data, ///< The new_cluster_data to query
		                                               const std::string      &arg_seq_name          ///< The name of the sequence of interest
		                                               ) {
			return has_domain_cluster_ids_of_seq_name( arg_new_cluster_data.get_dom_clust_ids(), arg_seq_name );
		}

		/// \brief Get the domain_cluster_ids associated with the specified sequence name in the specified new_cluster_data
		///
		/// \relates new_cluster_data
		inline const domain_cluster_ids & get_domain_cluster_ids_of_seq_name(const new_cluster_data &arg_new_cluster_data, ///< The new_cluster_data to query
		                                                                     const std::string      &arg_seq_name          ///< The name of the sequence of interest
		                                                                     ) {
			return get_domain_cluster_ids_of_seq_name( arg_new_cluster_data.get_dom_clust_ids(), arg_seq_name );
		}

		std::string to_string(const new_cluster_data &);

	} // namespace clust
} // namespace cath

#endif
