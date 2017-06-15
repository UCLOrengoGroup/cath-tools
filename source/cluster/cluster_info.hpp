/// \file
/// \brief The cluster_info class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_INFO_H
#define _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_INFO_H

#include <boost/utility/string_ref.hpp>

#include "cluster/cluster_name_ider.hpp"
#include "seq/seq_seg_run.hpp"

#include <string>

namespace cath {
	namespace clust {

		/// \brief Store a name/id lookup and size for clusters
		class cluster_info final {
		private:
			/// \brief The lookup between name and ID
			///
			/// The IDs for this can be used as to index into the sizes
			cluster_name_ider ider;

			/// \brief The sizes of the clusters
			size_vec sizes;

		public:
			/// \brief Default ctor
			cluster_info() = default;

			cluster_id_t add_name(const boost::string_ref &);
			const std::string & get_name_of_id(const size_t &);
			cluster_info & increment_size_of_cluster_of_id(const cluster_id_t &);
			size_t get_num_clusters() const;
			size_t get_size_of_cluster_of_id(const cluster_id_t &) const;
			const std::string & get_name_of_cluster_of_id(const cluster_id_t &) const;
		};

		/// \brief Add a cluster with the specified name
		///
		/// This can be called even if the name already exists
		///
		/// This sets the size for the cluster to 0. To add and increment in
		/// one call, use increment_and_get_id_for_cluster_of_name()
		inline cluster_id_t cluster_info::add_name(const boost::string_ref &arg_name ///< The name of the cluster to add
		                                           ) {
			const auto cluster_id = ider.add_name( std::string{ arg_name } );
			if ( cluster_id >= sizes.size() ) {
				sizes.resize( cluster_id + 1, 0_z );
			}
			return cluster_id;
		}

		/// \brief Get the name of the cluster with the specified ID
		inline const std::string & cluster_info::get_name_of_id(const size_t &arg_cluster_id ///< The ID of the cluster to get the name
		                                                        ) {
			return ider.get_name_of_id( arg_cluster_id );
		}

		/// \brief Increment the cluster with the specified ID
		inline cluster_info & cluster_info::increment_size_of_cluster_of_id(const cluster_id_t &arg_cluster_id ///< The ID of the cluster to increment the size of
		                                                                    ) {
			++( sizes[ arg_cluster_id ] );
			return *this;
		}

		/// \brief Get the number of clusters
		inline size_t cluster_info::get_num_clusters() const {
			return sizes.size();
		}

		/// \brief Get the size of the cluster with the specified ID
		inline size_t cluster_info::get_size_of_cluster_of_id(const cluster_id_t &arg_cluster_id ///< The ID of the cluster whose size should be retrieved
		                                                      ) const {
			return sizes[ arg_cluster_id ];
		}

		/// \brief Get the name of the cluster with the specified ID
		inline const std::string & cluster_info::get_name_of_cluster_of_id(const cluster_id_t &arg_cluster_id ///< The ID of the cluster whose name should be retrieved
		                                                                   ) const {
			return ider.get_name_of_id( arg_cluster_id );
		}

		/// \brief Increment the size for the cluster with the specified name in the specified cluster_infor
		///        and return its ID
		///
		/// \relates cluster_info
		inline cluster_id_t increment_and_get_id_for_cluster_of_name(cluster_info            &arg_cluster_info, ///< The cluster_info to query
		                                                             const boost::string_ref &arg_name          ///< The name of the cluster of interest
		                                                             ) {
			const auto cluster_id = arg_cluster_info.add_name( arg_name );
			arg_cluster_info.increment_size_of_cluster_of_id( cluster_id );
			return cluster_id;
		}

	} // namespace clust
} // namespace cath

#endif
