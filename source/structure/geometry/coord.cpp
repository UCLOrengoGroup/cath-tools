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

#include "coord.h"

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "exception/invalid_argument_exception.h"
#include "exception/out_of_range_exception.h"
#include "structure/geometry/angle.h"

using namespace cath::common;
using namespace cath::geom;
using namespace std;

using boost::lexical_cast;
using boost::property_tree::json_parser::write_json;
using boost::property_tree::ptree;

//const double LENGTH_CHECK_PRECISION_PERCENTAGE_TOLERANCE( 1E-10 );

const coord  coord::ORIGIN_COORD                        ( 0.0, 0.0, 0.0 );
const coord  coord::UNIT_X                              ( 1.0, 0.0, 0.0 );
const coord  coord::UNIT_Y                              ( 0.0, 1.0, 0.0 );
const coord  coord::UNIT_Z                              ( 0.0, 0.0, 1.0 );
const size_t coord::NUM_DIMS                            ( 3 );
const double coord::TOLERANCE_FOR_COORD_CLOSENESS_CHECKS( 0.00001 );

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

/// \brief TODOCUMENT
///
/// \relates coord
coord cath::geom::parallel_component_copy(const coord &arg_source_coord, ///< TODOCUMENT
                                          const coord &arg_dirn          ///< TODOCUMENT
                                          ) {
	const coord  unit_dirn = normalise_copy( arg_dirn );
	const double factor    = dot_product   ( arg_source_coord, unit_dirn );
	return factor * unit_dirn;
}

/// \brief TODOCUMENT
///
/// \relates coord
coord cath::geom::perpendicular_component_copy(const coord &arg_source_coord, ///< TODOCUMENT
                                               const coord &arg_dirn          ///< TODOCUMENT
                                               ) {
	const coord  unit_dirn = normalise_copy( arg_dirn );
	const double factor    = dot_product   ( arg_source_coord, unit_dirn );
	return arg_source_coord - ( factor * unit_dirn );
}

/// \brief TODOCUMENT
///
/// \relates coord
double cath::geom::length(const coord &arg_coord ///< TODOCUMENT
                          ) {
	return sqrt( squared_length( arg_coord ) );
}

/// \brief TODOCUMENT
///
/// \relates coord
double cath::geom::distance_between_points(const coord &arg_coord1, ///< TODOCUMENT
                                           const coord &arg_coord2  ///< TODOCUMENT
                                           ) {
	// Equivalently, either :
	//   length(arg_coord2 - arg_coord1)
	// or
	//   sqrt(squared_distance_between_points(arg_coord1, arg_coord2));
	return length(arg_coord2 - arg_coord1);
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
	return make_angle_from_radians<double>( acos( dot_product( normalise_copy( arg_vector_1 ), normalise_copy( arg_vector_2 ) ) ) );
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

/// \brief TODOCUMENT
///
/// This sends the
///
/// \relates coord
ostream & cath::geom::operator<<(ostream     &arg_os,   ///< TODOCUMENT
                                 const coord &arg_coord ///< TODOCUMENT
                                 ) {
//	ostringstream output_ss;
	arg_os << "coord[" << right << setw( 9 ) << arg_coord.get_x();
	arg_os << ", "     << right << setw( 9 ) << arg_coord.get_y();
	arg_os << ", "     << right << setw( 9 ) << arg_coord.get_z();
	arg_os << "]";
//	arg_os << output_ss.str();
	return arg_os;
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


/// \brief Make a new ptree representing the specified coord
///
/// \relates coord
ptree cath::geom::make_ptree_of(const coord &arg_coord ///< The coord that the new ptree should represent
                                ) {
	ptree new_ptree;
	save_to_ptree( new_ptree, arg_coord );
	return new_ptree;
}

/// \brief Create a JSON string to represent the specified coord
///
/// \relates coord
string cath::geom::to_json_string(const coord &arg_coord,       ///< The coord to represent in the JSON string
                                  const bool  &arg_pretty_print ///< Whether to use whitespace (including line breaks) in the JSON to make it more human-readable
                                  ) {
	ostringstream json_ss;
	ptree temp_ptree;
	save_to_ptree( temp_ptree, arg_coord );
	write_json( json_ss, temp_ptree, arg_pretty_print );
	return json_ss.str();
}
