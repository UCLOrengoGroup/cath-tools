/// \file
/// \brief The coord test suite

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

#include <array>
#include <string>
#include <utility>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/property_tree/from_json_string.hpp"
#include "cath/common/property_tree/to_json_string.hpp"
#include "cath/structure/geometry/angle.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/test/global_test_constants.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::geom;

using ::boost::lexical_cast;
using ::std::array;
using ::std::pair;
using ::std::string;

// clang-format off
constexpr coord coord1        = {  0.709591780486607204281312988314, -0.658089178396458751585385016369,  0.251789869421548129224674994475 };
constexpr coord coord2        = {  0.557086083496625361632936801470,  0.742781111328784104941291843716,  0.371391055664466840369186684256 };
constexpr coord coord3        = { -0.431433193716289575814215595528, -0.123267408225572736024666653520,  0.893683694284008178776446129632 };
constexpr coord double_coord1 = {  1.419183560973214408562625976630, -1.316178356792917503170770032740,  0.503579738843096258449349988950 };
// clang-format on

BOOST_FIXTURE_TEST_SUITE( coord_test_suite, cath::global_test_constants )

BOOST_AUTO_TEST_SUITE( string_conversion )

BOOST_AUTO_TEST_CASE( to_string_works ) {
	BOOST_CHECK_EQUAL( to_string( coord{ 800.0, 100.0, 400.0 } ), "coord[      800,       100,       400]" );
	BOOST_CHECK_EQUAL( to_string( coord{ 8.125, 1.125, 4.125 } ), "coord[    8.125,     1.125,     4.125]" );
}

BOOST_AUTO_TEST_CASE(lexical_cast_works) {
	BOOST_CHECK_EQUAL( lexical_cast<string>( coord{ 800.0, 100.0, 400.0 } ), "coord[      800,       100,       400]" );
	BOOST_CHECK_EQUAL( lexical_cast<string>( coord{ 8.125, 1.125, 4.125 } ), "coord[    8.125,     1.125,     4.125]" );
}

BOOST_AUTO_TEST_SUITE_END()

/// \brief Check equality
BOOST_AUTO_TEST_CASE( equality ) {
	static_assert( coord( 8.0, 1.0, 4.0 ) == coord( 8.0, 1.0, 4.0 ) );
	static_assert( coord( 5.0, 1.0, 4.0 ) != coord( 8.0, 1.0, 4.0 ) );
	static_assert( coord( 8.0, 5.0, 4.0 ) != coord( 8.0, 1.0, 4.0 ) );
	static_assert( coord( 8.0, 1.0, 5.0 ) != coord( 8.0, 1.0, 4.0 ) );

	BOOST_TEST( true );
}

/// \brief Check multiplication
BOOST_AUTO_TEST_CASE( multiplication ) {
	static_assert( double_coord1 == 2.0 * coord1 );
	static_assert( double_coord1 == coord1 * 2.0 );
	static_assert( double_coord1 == coord1 / 0.5 );
	static_assert( double_coord1 != 4.0 * coord1 );
	static_assert( double_coord1 != coord1 * 4.0 );
	static_assert( double_coord1 != coord1 / 0.25 );

	BOOST_TEST( true );
}

/// \brief Check lengths
BOOST_AUTO_TEST_CASE(lengths) {
	BOOST_CHECK_EQUAL(length(coord(0.0, 3.0, 4.0)), 5.0);
	BOOST_CHECK_EQUAL(length(coord(3.0, 4.0, 0.0)), 5.0);
	BOOST_CHECK_EQUAL(length(coord(4.0, 0.0, 3.0)), 5.0);
}

/// \brief Check components
BOOST_AUTO_TEST_CASE(components) {
	const coord coord_1_para_to_2 = parallel_component_copy(coord1, coord2);
	const coord coord_1_perp_to_2 = perpendicular_component_copy(coord1, coord2);
	BOOST_CHECK_EQUAL(coord1,         coord_1_para_to_2 + coord_1_perp_to_2);
	BOOST_CHECK_LT(                   dot_product(coord2,            coord_1_perp_to_2                ), ACCURACY_PERCENTAGE/100.0 );
	BOOST_CHECK_LT(                   dot_product(coord_1_para_to_2, coord_1_perp_to_2                ), ACCURACY_PERCENTAGE/100.0 );
	BOOST_CHECK_CLOSE(length(coord2), dot_product(coord2,            normalise_copy(coord_1_para_to_2)), ACCURACY_PERCENTAGE       );
}

/// \brief Check dot product of orthogonal vectors
BOOST_AUTO_TEST_CASE(dot_products) {
	BOOST_CHECK_CLOSE(dot_product(coord(1.41421356237, 1.41421356237, 0.0), coord(2.0, 0.0, 0.0)), 2.82842712475, ACCURACY_PERCENTAGE);
	static_assert(dot_product(coord(1.0, 0.0, 0.0), coord(1.0, 0.0, 0.0)) == 1.0);
	static_assert(dot_product(coord(0.0, 1.0, 0.0), coord(0.0, 1.0, 0.0)) == 1.0);
	static_assert(dot_product(coord(0.0, 0.0, 1.0), coord(0.0, 0.0, 1.0)) == 1.0);
	static_assert(dot_product(coord1, coord2) < ACCURACY_PERCENTAGE/100.0);
	static_assert(dot_product(coord2, coord3) < ACCURACY_PERCENTAGE/100.0);
	static_assert(dot_product(coord3, coord1) < ACCURACY_PERCENTAGE/100.0);
}

/// \brief Check cross_product
BOOST_AUTO_TEST_CASE(cross_products) {
	BOOST_CHECK_CLOSE(cross_product(coord1, coord2).get_x(), coord3.get_x(), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(cross_product(coord1, coord2).get_y(), coord3.get_y(), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(cross_product(coord1, coord2).get_z(), coord3.get_z(), ACCURACY_PERCENTAGE);

	BOOST_CHECK_CLOSE(cross_product(coord2, coord3).get_x(), coord1.get_x(), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(cross_product(coord2, coord3).get_y(), coord1.get_y(), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(cross_product(coord2, coord3).get_z(), coord1.get_z(), ACCURACY_PERCENTAGE);

	BOOST_CHECK_CLOSE(cross_product(coord3, coord1).get_x(), coord2.get_x(), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(cross_product(coord3, coord1).get_y(), coord2.get_y(), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(cross_product(coord3, coord1).get_z(), coord2.get_z(), ACCURACY_PERCENTAGE);

	BOOST_CHECK_CLOSE(length(cross_product(coord1, coord2)), length(coord3), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(length(cross_product(coord2, coord3)), length(coord1), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(length(cross_product(coord3, coord1)), length(coord2), ACCURACY_PERCENTAGE);

	BOOST_CHECK_CLOSE(length(cross_product(2.0 * coord1, 3.0 * coord2)), length(6.0 * coord3), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(length(cross_product(2.0 * coord2, 3.0 * coord3)), length(6.0 * coord1), ACCURACY_PERCENTAGE);
	BOOST_CHECK_CLOSE(length(cross_product(2.0 * coord3, 3.0 * coord1)), length(6.0 * coord2), ACCURACY_PERCENTAGE);
}
/// \brief Check angle_between_two_vectors() works as expected
BOOST_AUTO_TEST_CASE( angle_between_two_vectors_works ) {
	constexpr array coord_pairs = { pair{ coord( 5, 0, 0 ), coord( -5, 0, 0 ) },
		                            pair{ coord( 0, 5, 0 ), coord( 0, -5, 0 ) },
		                            pair{ coord( 0, 0, 5 ), coord( 0, 0, -5 ) } };
	/// \TODO Come constexpr-for, use here to make these tests constexpr
	for ( const auto &[ coord_1a, coord_1b ] : coord_pairs ) {
		BOOST_CHECK_EQUAL( ONE_REVOLUTION<double> / 2.0, angle_between_two_vectors( coord_1a, coord_1b ) );
		BOOST_CHECK_EQUAL( ONE_REVOLUTION<double> / 2.0, angle_between_two_vectors( coord_1b, coord_1a ) );

		for ( const auto &[ coord_2a, coord_2b ] : coord_pairs ) {
			if ( coord_1a != coord_2a ) {
				BOOST_CHECK_EQUAL( ONE_REVOLUTION<double> / 4.0, angle_between_two_vectors( coord_1a, coord_2a ) );
				BOOST_CHECK_EQUAL( ONE_REVOLUTION<double> / 4.0, angle_between_two_vectors( coord_1a, coord_2b ) );
				BOOST_CHECK_EQUAL( ONE_REVOLUTION<double> / 4.0, angle_between_two_vectors( coord_1b, coord_2a ) );
				BOOST_CHECK_EQUAL( ONE_REVOLUTION<double> / 4.0, angle_between_two_vectors( coord_1b, coord_2b ) );
			}
		}
	}
}

/// \brief Check dihedral_angle_between_four_points() works as expected
BOOST_AUTO_TEST_CASE(dihedral_angle_between_four_points_works) {
	BOOST_CHECK_CLOSE(
		angle_in_degrees(       - ONE_REVOLUTION<double> / 8.0 ),
		angle_in_degrees( dihedral_angle_between_four_points(coord(0, -1,  1), coord(0, 0, 0), coord(1, 0, 0), coord(1, 0, 1) ) ),
		ACCURACY_PERCENTAGE
	);
	BOOST_CHECK_CLOSE(
		angle_in_degrees(         ONE_REVOLUTION<double> / 8.0 ),
		angle_in_degrees( dihedral_angle_between_four_points(coord(0,  1,  1), coord(0, 0, 0), coord(1, 0, 0), coord(1, 0, 1) ) ),
		ACCURACY_PERCENTAGE
	);
	BOOST_CHECK_CLOSE(
		angle_in_degrees( - 3.0 * ONE_REVOLUTION<double> / 8.0 ),
		angle_in_degrees( dihedral_angle_between_four_points(coord(0, -1, -1), coord(0, 0, 0), coord(1, 0, 0), coord(1, 0, 1) ) ),
		ACCURACY_PERCENTAGE
	);
	BOOST_CHECK_CLOSE(
		angle_in_degrees(   3.0 * ONE_REVOLUTION<double> / 8.0 ),
		angle_in_degrees( dihedral_angle_between_four_points(coord(0,  1, -1), coord(0, 0, 0), coord(1, 0, 0), coord(1, 0, 1) ) ),
		ACCURACY_PERCENTAGE
	);
}

/// \brief Check that the dihedral angle gives the expected answer for a known tricky case
///
/// This example was flagged up by Natalie because one of her cath-ssaps was failing.
/// Trying to calculate the phi/psi angles of the PDB 1cn1, involves trying to calculate
/// the dihedral angle between these four atoms in residues 25 and 26 on chain A :
///
///     ATOM    184  N   ILE A  25      25.000  39.500   7.100  1.00 28.30
///     ATOM    185  CA  ILE A  25      24.100  38.300   7.300  1.00 28.30
///     ATOM    186  C   ILE A  25      24.900  37.100   7.500  1.00 28.30
///     ATOM    192  N   GLY A  26      24.300  35.900   7.700  1.00 18.50
///
/// The correct answer is ±180° but float rounding errors were leading to a dot-product value of:
///  -1.0000000000000002220446049250313080847263336181641
/// ...which was then passed to acos which then returned an angle of NaN.
///
/// Now, the dot-produt result gets clamped in [-1, 1] before the call to acos.
BOOST_AUTO_TEST_CASE(dihedral_angle_between_four_points_handles_float_difficulties) {
	BOOST_CHECK_EQUAL(
		angle_in_degrees( dihedral_angle_between_four_points(
			coord{ 25.0, 39.5, 7.1 },
			coord{ 24.1, 38.3, 7.3 },
			coord{ 24.9, 37.1, 7.5 },
			coord{ 24.3, 35.9, 7.7 }
		) ),
		-180.0
	);
}

BOOST_AUTO_TEST_CASE(planar_angles) {
	BOOST_CHECK_EQUAL( planar_angle_between( coord{ 0, 0, 1 }, coord{ 1, 0, 0 }, coord{ 0, 1, 0 } ), make_angle_from_degrees<double>( -90 ) );
	BOOST_CHECK_EQUAL( planar_angle_between( coord{ 0, 0, 1 }, coord{ 0, 1, 0 }, coord{ 1, 0, 0 } ), make_angle_from_degrees<double>(  90 ) );
}

BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_origin) {
	BOOST_CHECK_EQUAL( to_json_string( ORIGIN_COORD, json_style::COMPACT ), R"({"x":"0","y":"0","z":"0"})" "\n" );
}

BOOST_AUTO_TEST_CASE(to_json_string_works_for_unit_x) {
	BOOST_CHECK_EQUAL( to_json_string( UNIT_X_COORD,       json_style::COMPACT ), R"({"x":"1","y":"0","z":"0"})" "\n" );
}

BOOST_AUTO_TEST_CASE(to_json_string_works_for_unit_y) {
	BOOST_CHECK_EQUAL( to_json_string( UNIT_Y_COORD,       json_style::COMPACT ), R"({"x":"0","y":"1","z":"0"})" "\n" );
}

BOOST_AUTO_TEST_CASE(to_json_string_works_for_unit_z) {
	BOOST_CHECK_EQUAL( to_json_string( UNIT_Z_COORD,       json_style::COMPACT ), R"({"x":"0","y":"0","z":"1"})" "\n" );
}

BOOST_AUTO_TEST_CASE(from_json_string_works_for_origin) {
	BOOST_CHECK_EQUAL( from_json_string<coord>( R"({"x":"0","y":"0","z":"0"})" ), ORIGIN_COORD );
}

BOOST_AUTO_TEST_CASE(from_json_string_works_for_unit_x) {
	BOOST_CHECK_EQUAL( from_json_string<coord>( R"({"x":"1","y":"0","z":"0"})" ), UNIT_X_COORD       );
}

BOOST_AUTO_TEST_CASE(from_json_string_works_for_unit_y) {
	BOOST_CHECK_EQUAL( from_json_string<coord>( R"({"x":"0","y":"1","z":"0"})" ), UNIT_Y_COORD       );
}

BOOST_AUTO_TEST_CASE(from_json_string_works_for_unit_z) {
	BOOST_CHECK_EQUAL( from_json_string<coord>( R"({"x":"0","y":"0","z":"1"})" ), UNIT_Z_COORD       );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
