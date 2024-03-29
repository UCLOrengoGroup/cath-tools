/// \file
/// \brief The cluster_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTER_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTER_LIST_HPP

#include "cath/cluster/cluster_domains.hpp"
#include "cath/cluster/cluster_type_aliases.hpp"
#include "cath/common/container/id_of_string.hpp"
#include "cath/seq/seq_seg_run.hpp"

namespace cath::clust {

	/// \brief Store a list of cluster_domains entries, one for each cluster
	class cluster_list final {
	private:
		/// \brief The vector of cluster_domains entries, one for each cluster
		cluster_domains_vec cluster_seq_domains;

		/// \brief Ensure that there is a cluster_domains for the specified index and return it
		cluster_domains & ensure_and_get_cluster_domains_of_cluster_id(const size_t &prm_index ///< The index of the cluster_domains to ensure and return
		                                                               ) {
			if ( prm_index >= cluster_seq_domains.size() ) {
				cluster_seq_domains.resize( prm_index + 1 );
			}
			return cluster_seq_domains[ prm_index ];
		}

	public:
		/// \brief A const_iterator type alias as part of making this a range over cluster_domains entries
		using const_iterator = cluster_domains_vec_citr;

		/// \brief Add a domain with the specified (optional) segments in the specified cluster
		///        under the specified sequence ID
		clust_entry_problem add_domain_to_cluster(const cluster_id_t                  &prm_cluster_id, ///< The cluster ID of the domain to add
		                                          const common::id_of_string::id_type &prm_seq_id,     ///< The ID of the sequence on which the domain to add appears
		                                          seq::seq_seg_run_opt                 prm_segments    ///< The (optional) segments of the domain to add
		                                          ) {
			return ensure_and_get_cluster_domains_of_cluster_id( prm_cluster_id ).add_domain(
				prm_seq_id,
				::std::move( prm_segments ),
				prm_cluster_id
			);
		}

		/// \brief Return whether this is empty (ie stores info for no clusters)
		[[nodiscard]] bool empty() const {
			return cluster_seq_domains.empty();
		}

		/// \brief Return the number of cluster_domains entries (ie the number of clusters)
		[[nodiscard]] size_t size() const {
			return cluster_seq_domains.size();
		}

		/// \brief Get the cluster_domains associated with the sequence with the specified index
		const cluster_domains & operator[](const size_t &prm_index ///< The index of the cluster_domains to retrieve
		                                   ) const {
			return cluster_seq_domains[ prm_index ];
		}

		/// \brief Standard const begin() method, as part of making this into a range over the cluster_domains entries
		[[nodiscard]] const_iterator begin() const {
			return ::std::cbegin( cluster_seq_domains );
		}

		/// \brief Standard const end() method, as part of making this into a range over the cluster_domains entries
		[[nodiscard]] const_iterator end() const {
			return ::std::cend( cluster_seq_domains );
		}
	};

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTER_LIST_HPP
