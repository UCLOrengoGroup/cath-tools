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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_COORD_LINKAGE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_COORD_LINKAGE_HPP

namespace cath::geom {

	/// \brief Whether a coord can be used to link to further coord
	///        (eg when deciding which non-protein atoms to include near a domain)
	enum class coord_linkage : bool {
		ADD_AND_LINK, ///< The coord can be used to link to other coords
		ADD_ONLY      ///< The coord can only be added, not used to link to other coords
	};

} // namespace cath::geom

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_COORD_LINKAGE_HPP
