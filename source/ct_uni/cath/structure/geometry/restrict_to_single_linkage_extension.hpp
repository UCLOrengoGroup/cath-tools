/// \file
/// \brief The restrict_to_single_linkage_extension header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_GEOMETRY_RESTRICT_TO_SINGLE_LINKAGE_EXTENSION_HPP
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_GEOMETRY_RESTRICT_TO_SINGLE_LINKAGE_EXTENSION_HPP

#include "cath/structure/geometry/coord_linkage.hpp"
#include "cath/structure/structure_type_aliases.hpp"

namespace cath {
	namespace geom {

		void restrict_to_single_linkage_extension(coord_coord_linkage_pair_vec &,
		                                          const coord_coord_linkage_pair_vec_itr &,
		                                          const double &);

		coord_coord_linkage_pair_vec restrict_to_single_linkage_extension_copy(coord_coord_linkage_pair_vec,
		                                                                       const size_t &,
		                                                                       const double &);

	} // namespace geom
} // namespace cath

#endif
