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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_MAP_MAP_CLUSTERS_H
#define _CATH_TOOLS_SOURCE_CLUSTER_MAP_MAP_CLUSTERS_H

#include <boost/optional.hpp>

#include "cluster/cluster_type_aliases.hpp"
#include "common/type_aliases.hpp"

#include <iostream>
#include <string>

namespace cath { namespace clust { class clust_mapping_spec; } }
namespace cath { namespace clust { class new_cluster_data; } }
namespace cath { namespace clust { struct map_results; } }

namespace cath {
	namespace clust {

		map_results map_clusters(const old_cluster_data_opt &,
		                         const new_cluster_data &,
		                         const clust_mapping_spec &,
		                         const ostream_ref_opt & = boost::none );

		size_vec get_info_ordered_indices_of_unmapped_new_clusters(const potential_map_vec &,
		                                                           const new_cluster_data &);

	} // namespace clust
} // namespace cath

#endif
