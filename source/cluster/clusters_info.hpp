/// \file
/// \brief The clusters_info class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_CLUSTERS_INFO_HPP
#define _CATH_TOOLS_SOURCE_CLUSTER_CLUSTERS_INFO_HPP

#include <boost/utility/string_ref.hpp>

#include "cluster/cluster_info.hpp"
#include "cluster/cluster_type_aliases.hpp"
#include "common/container/id_of_str_bidirnl.hpp"
#include "seq/seq_seg_run.hpp"

#include <string>

namespace cath {
	namespace clust {

		/// \brief Store a name/id lookup and size for clusters
		class clusters_info final {
		private:
			/// \brief The lookup between name and ID
			///
			/// The IDs for this can be used as to index into the sizes
			common::id_of_str_bidirnl ider;

			/// \brief The info of the clusters
			cluster_info_vec infos;

		public:
			/// \brief Default ctor
			clusters_info() = default;

			cluster_id_t add_name(const boost::string_ref &);
			const std::string & get_name_of_id(const size_t &);
			clusters_info & update_info_for_cluster_of_id(const cluster_id_t &,
			                                              const boost::string_ref &,
			                                              const seq::seq_seg_run_opt &);
			size_t get_num_clusters() const;
			const cluster_info & get_info_of_cluster_of_id(const cluster_id_t &) const;
			const std::string & get_name_of_cluster_of_id(const cluster_id_t &) const;

			const common::id_of_str_bidirnl & get_ider() const;
		};

		/// \brief Add a cluster with the specified name
		///
		/// This can be called even if the name already exists
		///
		/// This sets the size for the cluster to 0. To add and increment in
		/// one call, use increment_and_get_id_for_cluster_of_name()
		inline cluster_id_t clusters_info::add_name(const boost::string_ref &prm_name ///< The name of the cluster to add
		                                            ) {
			const auto cluster_id = ider.add_name( prm_name );
			if ( cluster_id >= infos.size() ) {
				infos.resize( cluster_id + 1, cluster_info{} );
			}
			return cluster_id;
		}

		/// \brief Get the name of the cluster with the specified ID
		inline const std::string & clusters_info::get_name_of_id(const size_t &prm_cluster_id ///< The ID of the cluster to get the name
		                                                         ) {
			return ider.get_name_of_id( prm_cluster_id );
		}

		/// \brief Increment the cluster with the specified ID
		inline clusters_info & clusters_info::update_info_for_cluster_of_id(const cluster_id_t         &prm_cluster_id, ///< The ID of the cluster to increment the size of
		                                                                    const boost::string_ref    &prm_name,       ///< The name of the new entry in the cluster
		                                                                    const seq::seq_seg_run_opt &prm_segments    ///< The (optional) segments of the new entry in the cluster
		                                                                    ) {
			infos[ prm_cluster_id ].add_entry( prm_name, prm_segments );
			return *this;
		}

		/// \brief Get the number of clusters
		inline size_t clusters_info::get_num_clusters() const {
			return infos.size();
		}

		/// \brief Get the cluster_info of the cluster with the specified ID
		inline const cluster_info & clusters_info::get_info_of_cluster_of_id(const cluster_id_t &prm_cluster_id ///< The ID of the cluster whose size should be retrieved
		                                                                     ) const {
			return infos[ prm_cluster_id ];
		}

		/// \brief Get the name of the cluster with the specified ID
		inline const std::string & clusters_info::get_name_of_cluster_of_id(const cluster_id_t &prm_cluster_id ///< The ID of the cluster whose name should be retrieved
		                                                                    ) const {
			return ider.get_name_of_id( prm_cluster_id );
		}

		/// \brief Get the id_of_str_bidirnl used to ID the clusters
		inline const common::id_of_str_bidirnl & clusters_info::get_ider() const {
			return ider;
		}

		/// \brief Update the info for the specified cluster to account for a new domain with the specified name and
		///        (optional) segments in the specified clusters_info and return the cluster's ID
		///
		/// \relates clusters_info
		inline cluster_id_t update_info_and_get_id_for_cluster_of_name(clusters_info              &prm_clusters_info, ///< The clusters_info to query
		                                                               const boost::string_ref    &prm_cluster_name,  ///< The name of the cluster of interest
		                                                               const boost::string_ref    &prm_entry_name,    ///< The name of the entry of interest
		                                                               const seq::seq_seg_run_opt &prm_entry_segments ///< The (optional) segments of the entry to add
		                                                               ) {
			const auto cluster_id = prm_clusters_info.add_name( prm_cluster_name );
			prm_clusters_info.update_info_for_cluster_of_id( cluster_id, prm_entry_name, prm_entry_segments );
			return cluster_id;
		}

	} // namespace clust
} // namespace cath

#endif
