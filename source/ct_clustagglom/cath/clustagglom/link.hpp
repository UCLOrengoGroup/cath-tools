/// \file
/// \brief The link class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINK_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINK_HPP

#include "cath/clustagglom/clustagglom_type_aliases.hpp"

namespace cath::clust {

	/// \brief Represent a (half) link in a clustering to an item (from some other implied item)
	struct link final {
		/// \brief The index of the item linked to
		item_idx node;

		/// \brief The degree of dissimilarity between the two items
		strength dissim;

		link(const item_idx &,
		     const strength &);
	};

	/// \brief Construct from an index and a dissimilarity
	inline link::link(const item_idx &prm_node,  ///< The index of the item linked to
	                  const strength &prm_dissim ///< The degree of dissimilarity between the two items
	                  ) : node   { prm_node   },
	                      dissim { prm_dissim } {
	}

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINK_HPP
