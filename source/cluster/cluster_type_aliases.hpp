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

#include <vector>

namespace cath { namespace clust { struct domain_cluster_id; } }

namespace cath {
	namespace clust {

		/// \brief Type alias for the type used to index clusters
		using cluster_id_t = size_t;

		/// \brief Type alias for a vector of domain_cluster_id
		using domain_cluster_id_vec = std::vector<domain_cluster_id>;

	} // namespace clust
} // namespace cath

#endif
