/// \file
/// \brief The cross_covariance_matrix class definitions

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

#include "cross_covariance_matrix.hpp"

#include <boost/range/combine.hpp>

#include "cath/common/gsl/gsl_matrix_wrp.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/coord_list.hpp"

using namespace ::cath::geom::detail;

using ::boost::range::combine;

/// \brief Calculate the cross-covariance matrix for the specified coord_lists
///
/// For speeding up, could try using a std:array<double, 9> and a gsl_matrix_view of it
/// (but then don't use a subroutine for the array, else it'll get destroyed on return)
gsl_matrix_wrp cath::geom::detail::cross_covariance_matrix(const coord_list &prm_coords_a, ///< The first  list of coords to superpose
                                                           const coord_list &prm_coords_b  ///< The second list of coords to superpose
                                                           ) {
	gsl_matrix_wrp result{ 3, 3 };
	gsl_matrix_set_zero( result.get_ptr() );

	/// \TODO Come C++17 and structure bindings, use here
	for (const auto &coord_pair : combine( prm_coords_a, prm_coords_b ) ) {
		const coord &coord_a = coord_pair.get<0>();
		const coord &coord_b = coord_pair.get<1>();

		gsl_matrix_wrp_increment( result, 0, 0, coord_a.get_x() * coord_b.get_x() );
		gsl_matrix_wrp_increment( result, 0, 1, coord_a.get_x() * coord_b.get_y() );
		gsl_matrix_wrp_increment( result, 0, 2, coord_a.get_x() * coord_b.get_z() );

		gsl_matrix_wrp_increment( result, 1, 0, coord_a.get_y() * coord_b.get_x() );
		gsl_matrix_wrp_increment( result, 1, 1, coord_a.get_y() * coord_b.get_y() );
		gsl_matrix_wrp_increment( result, 1, 2, coord_a.get_y() * coord_b.get_z() );

		gsl_matrix_wrp_increment( result, 2, 0, coord_a.get_z() * coord_b.get_x() );
		gsl_matrix_wrp_increment( result, 2, 1, coord_a.get_z() * coord_b.get_y() );
		gsl_matrix_wrp_increment( result, 2, 2, coord_a.get_z() * coord_b.get_z() );
	}

	return result;
}
