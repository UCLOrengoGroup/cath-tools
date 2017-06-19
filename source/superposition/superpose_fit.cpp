/// \file
/// \brief The superpose_fit class definitions

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

#include "superpose_fit.hpp"

#include <boost/range/combine.hpp>

#include "common/gsl/get_determinant.hpp"
#include "common/gsl/gsl_matrix_wrp.hpp"
#include "common/gsl/gsl_vector_wrp.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/geometry/rotation.hpp"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace cath::common;
using namespace cath::geom;
using namespace cath::geom::detail;

using boost::range::combine;

/// \brief Calculate the cross-covariance matrix for the specified coord_lists
///
/// For speeding up, could try using a std:array<double, 9> and a gsl_matrix_view of it
/// (but then don't use a subroutine for the array, else it'll get destroyed on return)
gsl_matrix_wrp cross_covariance_matrix(const coord_list &arg_coords_a, ///< The first  list of coords to superpose
                                       const coord_list &arg_coords_b  ///< The second list of coords to superpose
                                       ) {
	gsl_matrix_wrp result{ 3, 3 };
	gsl_matrix_set_zero( result.get_ptr() );

	/// \TODO Come C++17 and structure bindings, use here
	for (const auto &coord_pair : combine( arg_coords_a, arg_coords_b ) ) {
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

/// \brief Find the rotation that, when applied to the first specified coord_list,
///        best superposes it onto the second specified coord list
///
/// \pre `arg_coords_a.size() == arg_coords_b.size()` else an invalid_argument_exception will be thrown
///
/// \pre Both arg_coords_a and arg_coords_b must be translated to have the their centres of gravity at the origin
///      else bad stuff might happen (most likely: meaningless results will be returned)
///
/// This uses the Kabsch algorithm, eg see https://en.wikipedia.org/wiki/Kabsch_algorithm
rotation cath::sup::superpose_fit_1st_to_2nd(const coord_list &arg_coords_a, ///< The first  list of coords to superpose onto the second
                                             const coord_list &arg_coords_b  ///< The second list of coords
                                             ) {
	// Check the sizes match
	if ( arg_coords_a.size() != arg_coords_b.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("This subroutine cannot fit lists of coordinates of different length"));
	}

	// Grab a cross-covariance matrix
	auto x_cov_mat = cross_covariance_matrix( arg_coords_a, arg_coords_b );

	// Do a singular value decomposition of the matrix into svd_left . S. svd_right ^T
	gsl_matrix_wrp  V    { 3, 3 };
	gsl_vector_wrp  S    { 3    };
	gsl_vector_wrp  work { 3    };
	gsl_vector_set_zero( work.get_ptr() );
	gsl_linalg_SV_decomp(
		x_cov_mat.get_ptr(),
		V.get_ptr(),
		S.get_ptr(),
		work.get_ptr()
	);
	gsl_matrix_wrp &svd_left  = x_cov_mat;
	gsl_matrix_wrp &svd_right = V;

	// Calculate svd_right . svd_left ^T into out
	gsl_matrix_wrp out { 3, 3 };
	gsl_blas_dgemm(
		CblasNoTrans,
		CblasTrans,
		1.0,
		svd_right.get_ptr(),
		svd_left.get_ptr(),
		0.0,
		out.get_ptr()
	);

	// Multiply the right column of svd_left by the determinant of out
	const double determinant = get_determinant( out );
	gsl_matrix_set( svd_left.get_ptr(), 0, 2, gsl_matrix_get( svd_left.get_ptr(), 0, 2 ) * determinant );
	gsl_matrix_set( svd_left.get_ptr(), 1, 2, gsl_matrix_get( svd_left.get_ptr(), 1, 2 ) * determinant );
	gsl_matrix_set( svd_left.get_ptr(), 2, 2, gsl_matrix_get( svd_left.get_ptr(), 2, 2 ) * determinant );

	// Recalculate svd_right . svd_left ^T into out
	gsl_blas_dgemm(
		CblasNoTrans,
		CblasTrans,
		1.0,
		svd_right.get_ptr(),
		svd_left.get_ptr(),
		0.0,
		out.get_ptr()
	);

	// Return a rotation built from out
	return rotation{
		gsl_matrix_get( out.get_ptr(), 0, 0 ),
		gsl_matrix_get( out.get_ptr(), 0, 1 ),
		gsl_matrix_get( out.get_ptr(), 0, 2 ),

		gsl_matrix_get( out.get_ptr(), 1, 0 ),
		gsl_matrix_get( out.get_ptr(), 1, 1 ),
		gsl_matrix_get( out.get_ptr(), 1, 2 ),

		gsl_matrix_get( out.get_ptr(), 2, 0 ),
		gsl_matrix_get( out.get_ptr(), 2, 1 ),
		gsl_matrix_get( out.get_ptr(), 2, 2 )
	};
}

/// \brief Find the rotation that, when applied to the second specified coord_list,
///        best superposes it onto the first specified coord list
///
/// \pre `arg_coords_a.size() == arg_coords_b.size()` else an invalid_argument_exception will be thrown
///
/// \pre Both arg_coords_a and arg_coords_b must be translated to have the their centres of gravity at the origin
///      else bad stuff might happen (most likely: meaningless results will be returned)
rotation cath::sup::superpose_fit_2nd_to_1st(const coord_list &arg_coords_a, ///< The first  list of coords
                                             const coord_list &arg_coords_b  ///< The second list of coords to superpose onto the first
                                             ) {
	return superpose_fit_1st_to_2nd( arg_coords_b, arg_coords_a );
}
