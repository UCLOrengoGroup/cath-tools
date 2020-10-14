/// \file
/// \brief The orient definitions

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

#include "orient.hpp"

// #include <boost/range/join.hpp>
#include <boost/range/numeric.hpp>

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/gsl/gsl_matrix_wrp.hpp"
#include "cath/common/gsl/gsl_vector_wrp.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/detail/cross_covariance_matrix.hpp"
#include "cath/structure/geometry/rotation.hpp"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::geom::detail;
using namespace ::cath::geom;

using ::boost::accumulate;
// using ::boost::range::join;

/// \brief Calculate the x and y coordinates of the later-weighted centre of gravity
///        (in which each coordinate is given a weight that increases linearly throughout the coords)
///
/// These orient functions aim to put the primary PCA direction on the x-axis and the secondary on the y-axis;
/// they then choose whether to apply further x, y or z flips by seeing which (if any) best pulls the later coords to
/// the positive x/y quadrant (which makes the result independent of the structures' initial orientations and more stable).
doub_doub_pair cath::geom::detail::x_and_y_of_later_weighted_cog(const coord_list &prm_coords ///< The coords being oriented
                                                                 ) {
	const double factor = 2.0 / ( static_cast<double>( prm_coords.size() ) * static_cast<double>( prm_coords.size() + 1 ) );
	const coord the_coord = accumulate(
		indices( prm_coords.size() ),
		ORIGIN_COORD,
		[&] (const coord &x, const size_t &y) {
			return x + ( static_cast<double>( y + 1_z ) * factor * prm_coords[ y ] );
		}
	);
	return { the_coord.get_x(), the_coord.get_y() };
}

/// \brief Get the orienting transformation for the specified (pre-superposed) coords
///
/// These orient functions aim to put the primary PCA direction on the x-axis and the secondary on the y-axis;
/// they then choose whether to apply further x, y or z flips by seeing which (if any) best pulls the later coords to
/// the positive x/y quadrant (which makes the result independent of the structures' initial orientations and more stable).
///
/// Note: the order of the coordinates matters
coord_rot_pair cath::geom::get_orienting_transformation(const coord_list &prm_core_coords ///< The coords (that have already been superposed) to orient
                                                        // const coord_list &prm_extra_coords ///< TODOCUMENT
                                                        ) {
	const auto all_coords = prm_core_coords;
	// const auto all_coords = join( prm_core_coords, prm_extra_coords );

	// Get the centre-of-gravity of the coords and calculate a centred bunch of coords
	const coord cog = centre_of_gravity( all_coords );
	const coord_list centred_all_coords = all_coords - cog;

	// Grab a cross-covariance matrix
	auto x_cov_mat = cross_covariance_matrix( centred_all_coords, centred_all_coords );

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

	// Use the results to generate an initial orientation rotation
	const rotation initial_orientation = rotation_to_x_axis_and_x_y_plane(
		coord {
			gsl_matrix_get( svd_left.get_ptr(), 0, 0 ),
			gsl_matrix_get( svd_left.get_ptr(), 1, 0 ),
			gsl_matrix_get( svd_left.get_ptr(), 2, 0 )
		},
		coord {
			gsl_matrix_get( svd_left.get_ptr(), 0, 1 ),
			gsl_matrix_get( svd_left.get_ptr(), 1, 1 ),
			gsl_matrix_get( svd_left.get_ptr(), 2, 1 )
		}
	);

	// Grab the x_and_y_of_later_weighted_cog() of the coords after the initial orientation is applied
	const auto x_y_of_wcog_of_init_orientn = x_and_y_of_later_weighted_cog(
		rotate_copy(
			initial_orientation,
			prm_core_coords - cog
		)
	);

	// Build a flip-rotation to apply to x_and_y_of_later_weighted_cog so its x_and_y_of_later_weighted_cog()
	// is positive in both dimensions
	const rotation flip_rotn{
		copysign( 1.0, x_y_of_wcog_of_init_orientn.first ), 0.0, 0.0,
		0.0, copysign( 1.0, x_y_of_wcog_of_init_orientn.second ), 0.0,
		0.0, 0.0, copysign( 1.0, x_y_of_wcog_of_init_orientn.first * x_y_of_wcog_of_init_orientn.second )
	};

	// Return the centring translation and the (appropriately flipped) orientation rotation
	return {
		-cog,
		flip_rotn * initial_orientation
	};

}
