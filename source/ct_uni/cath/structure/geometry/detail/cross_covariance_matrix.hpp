/// \file
/// \brief The cross_covariance_matrix header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_DETAIL_CROSS_COVARIANCE_MATRIX_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_DETAIL_CROSS_COVARIANCE_MATRIX_HPP

// clang-format off
namespace cath::geom { class coord_list; }
namespace cath::geom::detail { class gsl_matrix_wrp; }
// clang-format on

namespace cath::geom::detail {

	geom::detail::gsl_matrix_wrp cross_covariance_matrix(const geom::coord_list &,
	                                                     const geom::coord_list &);

} // namespace cath::geom::detail
#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_DETAIL_CROSS_COVARIANCE_MATRIX_HPP
