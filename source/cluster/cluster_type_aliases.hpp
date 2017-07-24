/// \file
/// \brief The cluster type_aliases header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_TYPE_ALIASES_H
#define _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_TYPE_ALIASES_H

#include <boost/optional/optional_fwd.hpp>

#include <cstddef>
#include <vector>

namespace cath { namespace clust { class cluster_domains; } }
namespace cath { namespace clust { class cluster_info; } }
namespace cath { namespace clust { class cluster_list; } }
namespace cath { namespace clust { class old_cluster_data; } }
namespace cath { namespace clust { namespace detail { class mapping_job; } } }
namespace cath { namespace clust { namespace detail { struct seq_id_and_domain_cluster_ids_pair; } } }
namespace cath { namespace clust { struct domain_cluster_id; } }

namespace cath {
	namespace clust {
		
		/// \brief Type alias for a vector of cluster_domains entries
		using cluster_domains_vec      = std::vector<cluster_domains>;

		/// \brief Type alias for cluster_domains_vec's const_iterator type
		using cluster_domains_vec_citr = cluster_domains_vec::const_iterator;

		/// \brief Type alias for the type used to index clusters
		using cluster_id_t             = size_t;

		/// \brief Type alias for a vector of cluster_info objects
		using cluster_info_vec         = std::vector<cluster_info>;

		/// \brief Type alias for a vector of domain_cluster_id
		using domain_cluster_id_vec    = std::vector<domain_cluster_id>;

		/// \brief Type alias for an optional old_cluster_data
		using old_cluster_data_opt     = boost::optional<old_cluster_data>;

		namespace detail {

			/// \brief Type alias for a vector of mapping_job objects
			using mapping_job_vec = std::vector<mapping_job>;

			/// \brief Type alias for a vector of seq_id_and_domain_cluster_ids_pair objects
			using seq_id_and_domain_cluster_ids_pair_vec      = std::vector<seq_id_and_domain_cluster_ids_pair>;

			/// \brief Type alias for seq_id_and_domain_cluster_ids_pair_vec's const_iterator type
			using seq_id_and_domain_cluster_ids_pair_vec_citr = seq_id_and_domain_cluster_ids_pair_vec::const_iterator;

		} // namespace detail

	} // namespace clust
} // namespace cath

#endif
