/// \file
/// \brief The rotation test suite

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
/// You should have received A copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <boost/lexical_cast.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/property_tree/from_json_string.hpp"
#include "cath/common/property_tree/to_json_string.hpp"
#include "cath/structure/geometry/angle.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/rotation.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "cath/test/global_test_constants.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::geom;

constexpr coord A = { 3, 7, 9 };
constexpr coord B = { 8, 5, 1 };

BOOST_FIXTURE_TEST_SUITE(rotation_test_suite, global_test_constants)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(rotation_to_x_axis_and_x_y_plane_works) {
	const rotation test_rotation = rotation_to_x_axis_and_x_y_plane(A, B);
	const coord    rotated_a     = rotate_copy(test_rotation, A);
	const coord    rotated_b     = rotate_copy(test_rotation, B);

	// Check that A's length has been preserved and that it has been rotated onto the x axis
	BOOST_CHECK_CLOSE( length(A), length(rotated_a),                             100.0 * rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS );
	BOOST_CHECK_CLOSE( length(A), dot_product( coord(1.0, 0.0, 0.0), rotated_a), 100.0 * rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS );

	// Check that B's length has been preserved and that it has been rotated onto the x-y plane
	BOOST_CHECK_EQUAL( length(B), length(rotated_b) );
	BOOST_CHECK_LE   ( dot_product( coord(0.0, 0.0, 1.0), rotated_b), rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS );
}

BOOST_AUTO_TEST_CASE(rotate_x_to_y_to_z_to_x_is_correct) {
	BOOST_TEST( rotate_copy( ROTATE_X_TO_Y_TO_Z_TO_X, UNIT_X_COORD ) == UNIT_Y_COORD );
	BOOST_TEST( rotate_copy( ROTATE_X_TO_Y_TO_Z_TO_X, UNIT_Y_COORD ) == UNIT_Z_COORD );
	BOOST_TEST( rotate_copy( ROTATE_X_TO_Y_TO_Z_TO_X, UNIT_Z_COORD ) == UNIT_X_COORD );
}

BOOST_AUTO_TEST_CASE(rotate_x_to_z_to_y_to_x_is_correct) {
	BOOST_TEST( rotate_copy( ROTATE_X_TO_Z_TO_Y_TO_X, UNIT_X_COORD ) == UNIT_Z_COORD );
	BOOST_TEST( rotate_copy( ROTATE_X_TO_Z_TO_Y_TO_X, UNIT_Z_COORD ) == UNIT_Y_COORD );
	BOOST_TEST( rotate_copy( ROTATE_X_TO_Z_TO_Y_TO_X, UNIT_Y_COORD ) == UNIT_X_COORD );
}

BOOST_AUTO_TEST_CASE(identity_rotation_is_correct) {
	BOOST_TEST( rotate_copy( IDENTITY_ROTATION, UNIT_X_COORD ) == UNIT_X_COORD );
	BOOST_TEST( rotate_copy( IDENTITY_ROTATION, UNIT_Y_COORD ) == UNIT_Y_COORD );
	BOOST_TEST( rotate_copy( IDENTITY_ROTATION, UNIT_Z_COORD ) == UNIT_Z_COORD );
}

/// \brief Check that this example (1c0pA01, residue 999) can be used to produce A similar but more accurate rotation
BOOST_AUTO_TEST_CASE(tidy_rotation_does_not_always_throw) {
	BOOST_CHECK_NO_THROW_DIAG(
		tidy_rotation( 0.0631, -0.9000, -0.4312,
		               0.9437,  0.1944, -0.2676,
		               0.3247, -0.3901,  0.8616,
		               0.0001
		)
	);
}

/// \brief Check that an attempt to find A similar but more rotation for this example (1c0pA01, residue 999) can
///        fail if the tolerance is too strict
BOOST_AUTO_TEST_CASE(tidy_rotation_can_fail) {
	BOOST_CHECK_THROW(
		tidy_rotation( 0.0631, -0.9000, -0.4312,
		               0.9437,  0.1944, -0.2676,
		               0.3247, -0.3901,  0.8616,
		               0.000001
		),
		invalid_argument_exception
	);
}

/// \brief Check that tidy_copy() returns exact copies of standard, tidy rotations
BOOST_AUTO_TEST_CASE(tidy_copy_works) {
	BOOST_CHECK_EQUAL( IDENTITY_ROTATION,       tidy_copy( IDENTITY_ROTATION       ) );
	BOOST_CHECK_EQUAL( ROTATE_X_TO_Y_TO_Z_TO_X, tidy_copy( ROTATE_X_TO_Y_TO_Z_TO_X ) );
	BOOST_CHECK_EQUAL( ROTATE_X_TO_Z_TO_Y_TO_X, tidy_copy( ROTATE_X_TO_Z_TO_Y_TO_X ) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(rotation_angle_works) {
	BOOST_CHECK_EQUAL(                               ZERO_ANGLE<double>,                   angle_of_rotation( IDENTITY_ROTATION       )                        );
	BOOST_CHECK_CLOSE( angle_in_degrees( ONE_REVOLUTION<double> / 3.0 ), angle_in_degrees( angle_of_rotation( ROTATE_X_TO_Y_TO_Z_TO_X ) ), ACCURACY_PERCENTAGE );
	BOOST_CHECK_CLOSE( angle_in_degrees( ONE_REVOLUTION<double> / 3.0 ), angle_in_degrees( angle_of_rotation( ROTATE_X_TO_Z_TO_Y_TO_X ) ), ACCURACY_PERCENTAGE );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(angle_between_rotations_works) {
	BOOST_CHECK_EQUAL(                               ZERO_ANGLE<double>,                   angle_between_rotations( IDENTITY_ROTATION,       IDENTITY_ROTATION          )                        );
	BOOST_CHECK_CLOSE( angle_in_degrees( ONE_REVOLUTION<double> / 3.0 ), angle_in_degrees( angle_between_rotations( IDENTITY_ROTATION,       ROTATE_X_TO_Y_TO_Z_TO_X    ) ), ACCURACY_PERCENTAGE );
	BOOST_CHECK_CLOSE( angle_in_degrees( ONE_REVOLUTION<double> / 3.0 ), angle_in_degrees( angle_between_rotations( IDENTITY_ROTATION,       ROTATE_X_TO_Z_TO_Y_TO_X    ) ), ACCURACY_PERCENTAGE );

	BOOST_CHECK_CLOSE( angle_in_degrees( ONE_REVOLUTION<double> / 3.0 ), angle_in_degrees( angle_between_rotations( ROTATE_X_TO_Y_TO_Z_TO_X, IDENTITY_ROTATION          ) ), ACCURACY_PERCENTAGE );
	BOOST_CHECK_EQUAL(                               ZERO_ANGLE<double>,                   angle_between_rotations( ROTATE_X_TO_Y_TO_Z_TO_X, ROTATE_X_TO_Y_TO_Z_TO_X    )                        );
	BOOST_CHECK_CLOSE( angle_in_degrees( ONE_REVOLUTION<double> / 3.0 ), angle_in_degrees( angle_between_rotations( ROTATE_X_TO_Y_TO_Z_TO_X, ROTATE_X_TO_Z_TO_Y_TO_X    ) ), ACCURACY_PERCENTAGE );

	BOOST_CHECK_CLOSE( angle_in_degrees( ONE_REVOLUTION<double> / 3.0 ), angle_in_degrees( angle_between_rotations( ROTATE_X_TO_Z_TO_Y_TO_X, IDENTITY_ROTATION          ) ), ACCURACY_PERCENTAGE );
	BOOST_CHECK_CLOSE( angle_in_degrees( ONE_REVOLUTION<double> / 3.0 ), angle_in_degrees( angle_between_rotations( ROTATE_X_TO_Z_TO_Y_TO_X, ROTATE_X_TO_Y_TO_Z_TO_X    ) ), ACCURACY_PERCENTAGE );
	BOOST_CHECK_EQUAL(                               ZERO_ANGLE<double>,                   angle_between_rotations( ROTATE_X_TO_Z_TO_Y_TO_X, ROTATE_X_TO_Z_TO_Y_TO_X    )                        );

	BOOST_CHECK_EQUAL(                   ONE_REVOLUTION<double> / 2.0  ,                   angle_between_rotations( IDENTITY_ROTATION,       rotation(0, -1, 0, -1, 0, 0, 0, 0, -1) )                          );
	BOOST_CHECK_EQUAL(                   ONE_REVOLUTION<double> / 2.0  ,                   angle_between_rotations( ROTATE_X_TO_Y_TO_Z_TO_X, rotation(0, -1, 0, -1, 0, 0, 0, 0, -1) )                          );
	BOOST_CHECK_EQUAL(                   ONE_REVOLUTION<double> / 2.0  ,                   angle_between_rotations( ROTATE_X_TO_Z_TO_Y_TO_X, rotation(0, -1, 0, -1, 0, 0, 0, 0, -1) )                          );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(rotation_of_angle_works) {
	const doub_angle_vec angles = { make_angle_from_revolutions<double>( 0.000 ),
	                                make_angle_from_revolutions<double>( 0.100 ),
	                                make_angle_from_revolutions<double>( 0.250 ),
	                                make_angle_from_revolutions<double>( 0.333 ),
	                                make_angle_from_revolutions<double>( 0.500 ),
	                                make_angle_from_degrees<double>    ( 1.000 ),
	                                make_angle_from_degrees<double>    ( 2.000 ),
	                                make_angle_from_radians<double>    ( 1.000 ),
	                                make_angle_from_radians<double>    ( 2.000 ) };
	for (const auto &the_angle : angles) {
		BOOST_CHECK_CLOSE( angle_in_degrees( the_angle ), angle_in_degrees( angle_of_rotation( rotation_of_angle( the_angle ) ) ), ACCURACY_PERCENTAGE );
	}
}

BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_SUITE(write)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_identity) {
	BOOST_CHECK_EQUAL( to_json_string( IDENTITY_ROTATION,       json_style::COMPACT ), R"({"":["1","0","0"],"":["0","1","0"],"":["0","0","1"]})" "\n" );
}

BOOST_AUTO_TEST_CASE(to_json_string_works_for_x_to_y_to_z_to_x) {
	BOOST_CHECK_EQUAL( to_json_string( ROTATE_X_TO_Y_TO_Z_TO_X, json_style::COMPACT ), R"({"":["0","0","1"],"":["1","0","0"],"":["0","1","0"]})" "\n" );
}

BOOST_AUTO_TEST_CASE(to_json_string_works_for_x_to_z_to_y_to_x) {
	BOOST_CHECK_EQUAL( to_json_string( ROTATE_X_TO_Z_TO_Y_TO_X, json_style::COMPACT ), R"({"":["0","1","0"],"":["0","0","1"],"":["1","0","0"]})" "\n" );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(read)

BOOST_AUTO_TEST_CASE(throws_on_from_json_string_with_missing_row) {
	BOOST_CHECK_THROW( from_json_string<rotation>( R"({"":["1","0","0"],"":["0","1","0"],"":["0","0"]})"          "\n" ), runtime_error_exception );
}

BOOST_AUTO_TEST_CASE(throws_on_from_json_string_with_short_row) {
	BOOST_CHECK_THROW( from_json_string<rotation>( R"({"":["1","0","0"],"":["0","1","0"],"":["0","0"]})"          "\n" ), runtime_error_exception );
}

BOOST_AUTO_TEST_CASE(throws_on_from_json_string_with_spurious_label) {
	BOOST_CHECK_THROW( from_json_string<rotation>( R"({"":["1","0","0"],"":["0","1","0"],"wrong":["0","0","1"]})" "\n" ), runtime_error_exception );
}

BOOST_AUTO_TEST_CASE(from_json_string_works_for_identity) {
	BOOST_CHECK_EQUAL( from_json_string<rotation>( R"({"":["1","0","0"],"":["0","1","0"],"":["0","0","1"]})" "\n" ), IDENTITY_ROTATION       );
}

BOOST_AUTO_TEST_CASE(from_json_string_works_for_x_to_y_to_z_to_x) {
	BOOST_CHECK_EQUAL( from_json_string<rotation>( R"({"":["0","0","1"],"":["1","0","0"],"":["0","1","0"]})" "\n" ), ROTATE_X_TO_Y_TO_Z_TO_X );
}

BOOST_AUTO_TEST_CASE(from_json_string_works_for_x_to_z_to_y_to_x) {
	BOOST_CHECK_EQUAL( from_json_string<rotation>( R"({"":["0","1","0"],"":["0","0","1"],"":["1","0","0"]})" "\n" ), ROTATE_X_TO_Z_TO_Y_TO_X );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

