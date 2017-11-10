/// \file
/// \brief The coord class definitions

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

#include "coord.hpp"

#include <boost/algorithm/clamp.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "structure/geometry/angle.hpp"

using namespace cath::common;
using namespace cath::geom;
using namespace std;

using boost::algorithm::clamp;
using boost::format;
using boost::property_tree::ptree;

//const double LENGTH_CHECK_PRECISION_PERCENTAGE_TOLERANCE( 1E-10 );
constexpr coord  coord::ORIGIN_COORD { 0.0, 0.0, 0.0 };
constexpr coord  coord::UNIT_X       { 1.0, 0.0, 0.0 };
constexpr coord  coord::UNIT_Y       { 0.0, 1.0, 0.0 };
constexpr coord  coord::UNIT_Z       { 0.0, 0.0, 1.0 };
constexpr size_t coord::NUM_DIMS;
constexpr double coord::TOLERANCE_FOR_COORD_CLOSENESS_CHECKS;

/// \brief TODOCUMENT
///
/// \relates coord
coord cath::geom::normalise_copy(const coord &arg_coord ///< TODOCUMENT
                                 ) {
	const double coord_length( length( arg_coord ) );

#ifndef NDEBUG
	using boost::math::isnormal;
	if ( ! isnormal( coord_length ) ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Invalid coord factor"));
	}
#endif
	return ( arg_coord / coord_length );
}


/// \brief Compute the angle (in radians) between two vectors (represented as coord objects)
///
/// \relates coord
///
/// If either vector is of zero length, an invalid_argument_exception will be thrown.
doub_angle cath::geom::angle_between_two_vectors(const coord &arg_vector_1, ///< TODOCUMENT
                                                 const coord &arg_vector_2  ///< TODOCUMENT
                                                 ) {
#ifndef NDEBUG
	using boost::math::isnormal;
	if ( ! isnormal( length( arg_vector_1 ) ) || ! isnormal( length( arg_vector_2 ) ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("AngleBetweenZeroLengthVectorsRequested"));
	}
#endif
	return make_angle_from_radians<double>( acos(
		clamp(
			dot_product(
				normalise_copy( arg_vector_1 ),
				normalise_copy( arg_vector_2 )
			),
			-1.0,
			 1.0
		)
	) );
}

/// \brief Compute the angle (in radians) determined by three points
///
/// \relates coord
///
/// Three points must be given as arguments and the result is the angle between the two straight lines joining 1 to 2 and 2 to 3.
/// If point 2 is the same as either point 1 or point 2 (or both), an invalid_argument_exception will be thrown.
doub_angle cath::geom::angle_between_three_points(const coord &arg_coord1, ///< TODOCUMENT
                                                  const coord &arg_coord2, ///< TODOCUMENT
                                                  const coord &arg_coord3  ///< TODOCUMENT
                                                  ) {
	const coord vector_to_1( arg_coord1 - arg_coord2 );
	const coord vector_to_3( arg_coord3 - arg_coord2 );
	return angle_between_two_vectors( vector_to_1, vector_to_3 );
}

/// \brief TODOCUMENT
///
/// \relates coord
doub_angle cath::geom::dihedral_angle_between_four_points(const coord &arg_coord1, ///< TODOCUMENT
                                                          const coord &arg_coord2, ///< TODOCUMENT
                                                          const coord &arg_coord3, ///< TODOCUMENT
                                                          const coord &arg_coord4  ///< TODOCUMENT
                                                          ) {
	const coord vector_1_to_2( arg_coord2 - arg_coord1 );
	const coord vector_2_to_3( arg_coord3 - arg_coord2 );
	const coord vector_3_to_4( arg_coord4 - arg_coord3 );

	const coord unit_perp_vect_at_2( normalise_copy( cross_product( vector_1_to_2, vector_2_to_3 ) ) );
	const coord unit_perp_vect_at_3( normalise_copy( cross_product( vector_2_to_3, vector_3_to_4 ) ) );

	const auto   the_angle = angle_between_two_vectors( unit_perp_vect_at_2, unit_perp_vect_at_3 );
	const double sign      = dot_product( vector_2_to_3, cross_product( unit_perp_vect_at_2, unit_perp_vect_at_3 ) );

	return ( ( sign >= 0) ? the_angle : -the_angle );
}

/// \brief Calculate the angle between two specified coords when projected onto a plane orthogoanl to a specified coord
///
/// When looking at the plane along the direction of arg_plane_ortho, the angle
/// is from the projection of arg_coord1 to the projection of arg_coord2,
/// where (like on the complex plane) anti-clockwise is positive.
///
/// This is used to produces sec files
///
/// \relates coord
doub_angle cath::geom::planar_angle_between(const coord &arg_plane_ortho, ///< The coord defining the plane
                                            const coord &arg_coord1,      ///< The first  angle to compare
                                            const coord &arg_coord2       ///< The second angle to compare
                                            ) {
	const coord      in_plane1 = normalise_copy( perpendicular_component_copy( arg_coord1, arg_plane_ortho ) );
	const coord      in_plane2 = normalise_copy( perpendicular_component_copy( arg_coord2, arg_plane_ortho ) );
	const double     sign      = dot_product( cross_product( in_plane2, in_plane1 ), normalise_copy( arg_plane_ortho ) );
	const doub_angle ang_rads  = angle_between_two_vectors(
		perpendicular_component_copy( arg_coord1, arg_plane_ortho ),
		perpendicular_component_copy( arg_coord2, arg_plane_ortho )
	);
	return ( sign >= 0.0 ) ? ang_rads : -ang_rads;
}

/// \brief Generate a string describing the specified coord
///
/// \relates coord
string cath::geom::to_string(const coord &arg_coord ///< The coord to describe
                             ) {
	return "coord["
		+ ( format( R"(%9g)"  ) % arg_coord.get_x() ).str()
		+ ", "
		+ ( format( R"(%9g)" ) % arg_coord.get_y() ).str()
		+ ", "
		+ ( format( R"(%9g)"  ) % arg_coord.get_z() ).str()
		+ "]";
}

/// \brief Insert a description of the specified coord into the specified ostream
///
/// \relates coord
ostream & cath::geom::operator<<(ostream     &arg_os,   ///< The ostream into which the description should be inserted
                                 const coord &arg_coord ///< The coord to describe
                                 ) {
	arg_os << to_string( arg_coord );
	return arg_os;
}

/// \brief Build a coord from a coord-populated ptree
///
/// \relates coord
coord cath::geom::coord_from_ptree(const ptree &arg_ptree ///< The ptree from which the coord should be read
                                   ) {
	return {
		arg_ptree.get<double>( "x" ),
		arg_ptree.get<double>( "y" ),
		arg_ptree.get<double>( "z" )
	};
}

/// \brief Save the specified coord to the specified Boost Property Tree ptree
///
/// \relates coord
void cath::geom::save_to_ptree(ptree       &arg_ptree, ///< The ptree to which the coord should be saved
                               const coord &arg_coord  ///< The coord to be saved
                               ) {
	arg_ptree.put( "x", arg_coord.get_x() );
	arg_ptree.put( "y", arg_coord.get_y() );
	arg_ptree.put( "z", arg_coord.get_z() );
}
