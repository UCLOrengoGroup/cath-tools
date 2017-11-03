/// \file
/// \brief The coord_linkage header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_COORD_LINKAGE_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_COORD_LINKAGE_H

#include <utility>
#include <vector>

namespace cath { namespace geom { class coord; } }

namespace cath {
	namespace geom {

		/// \brief Whether a coord can be used to link to further coord
		///        (eg when deciding which non-protein atoms to include near a domain)
		enum class coord_linkage : bool {
			ADD_AND_LINK, ///< The coord can be used to link to other coords
			ADD_ONLY      ///< The coord can only be added, not used to link to other coords
		};

		/// \brief Type alias for a pair of coord and coord_linkage
		///
		/// \todo Consider moving into structure_type_aliases.hpp
		/// (and then remove the #include statements from this header)
		using coord_coord_linkage_pair = std::pair<coord, coord_linkage>;

		/// \brief Type alias for a vector of coord_coord_linkage_pair values
		///
		/// \todo Consider moving into structure_type_aliases.hpp
		/// (and then remove the #include statements from this header)
		using coord_coord_linkage_pair_vec = std::vector<coord_coord_linkage_pair>;

		/// \brief Type alias for coord_coord_linkage_pair_vec's iterator type
		///
		/// \todo Consider moving into structure_type_aliases.hpp
		/// (and then remove the #include statements from this header)
		using coord_coord_linkage_pair_vec_itr = coord_coord_linkage_pair_vec::iterator;

	} // namespace geom
} // namespace cath

#endif
