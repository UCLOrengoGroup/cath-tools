/// \file
/// \brief The eigen class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_PCA_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_PCA_H

#include "common/type_aliases.hpp"

namespace cath { namespace geom { class coord; } }
namespace cath { namespace geom { class coord_list; } }
namespace cath { namespace geom { class line; } }
namespace cath { namespace geom { namespace detail { class gsl_matrix_wrp; } } }
namespace cath { namespace geom { using coord_coord_pair = std::pair<coord, coord>; } } ///< \todo Move into structure/structure_type_aliases.hpp

namespace cath {
	namespace geom {

		namespace detail {
			doub_vec build_matrix_of_coords(const coord_list &,
			                                const coord &);
		} // namespace detail

		line line_of_best_fit(const coord_list &);

	} // namespace geom
} // namespace cath

#endif
