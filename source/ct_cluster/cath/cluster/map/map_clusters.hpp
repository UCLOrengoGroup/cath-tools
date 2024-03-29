/// \map
/// \brief The map_clusters class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_MAP_MAP_CLUSTERS_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_MAP_MAP_CLUSTERS_HPP

#include <iostream>
#include <optional>
#include <string>

#include "cath/cluster/cluster_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::clust { class clust_mapping_spec; }
namespace cath::clust { class new_cluster_data; }
namespace cath::clust { struct map_results; }
// clang-format on

namespace cath::clust {

	map_results map_clusters(const old_cluster_data_opt &,
	                         const new_cluster_data &,
	                         const clust_mapping_spec &,
	                         const ostream_ref_opt & = ::std::nullopt );

	size_vec get_info_ordered_indices_of_unmapped_new_clusters(const potential_map_vec &,
	                                                           const new_cluster_data &);

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_MAP_MAP_CLUSTERS_HPP
