/// \file
/// \brief The orient header

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
/// MERCHANTABILITY or ORIENTNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ORIENT_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ORIENT_HPP

#include "cath/common/type_aliases.hpp"
#include "cath/structure/structure_type_aliases.hpp"

namespace cath { namespace geom { class coord_list; } }
namespace cath { namespace geom { namespace detail { class gsl_matrix_wrp; } } }

namespace cath {
	namespace geom {
		namespace detail {

			doub_doub_pair x_and_y_of_later_weighted_cog(const coord_list &a);

		} // namespace detail

		coord_rot_pair get_orienting_transformation(const coord_list &);

	} // namespace geom
} // namespace cath
#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ORIENT_HPP
