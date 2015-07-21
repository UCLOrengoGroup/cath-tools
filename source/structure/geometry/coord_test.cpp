/// \file
/// \brief The coord test suite

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include <boost/math/constants/constants.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "structure/geometry/angle.h"
#include "structure/geometry/coord.h"
#include "test/global_test_constants.h"

using namespace boost::math::constants;
using namespace cath;
using namespace cath::geom;
using namespace std;

namespace cath {
	namespace test {

		struct coord_test_suite_fixture : public global_test_constants {
		protected:
			~coord_test_suite_fixture() noexcept = default;

			const coord coord1        = {  0.709591780486607204281312988314, -0.658089178396458751585385016369,  0.251789869421548129224674994475 };
			const coord coord2        = {  0.557086083496625361632936801470,  0.742781111328784104941291843716,  0.371391055664466840369186684256 };
			const coord coord3        = { -0.431433193716289575814215595528, -0.123267408225572736024666653520,  0.893683694284008178776446129632 };
			const coord double_coord1 = {  1.419183560973214408562625976630, -1.316178356792917503170770032740,  0.503579738843096258449349988950 };
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(coord_test_suite, cath::test::coord_test_suite_fixture)

/// \brief Check equality
BOOST_AUTO_TEST_CASE(equality) {
	BOOST_CHECK_EQUAL( coord(8.0, 1.0, 4.0), coord(8.0, 1.0, 4.0) );
	BOOST_CHECK_NE(    coord(5.0, 1.0, 4.0), coord(8.0, 1.0, 4.0) );
	BOOST_CHECK_NE(    coord(8.0, 5.0, 4.0), coord(8.0, 1.0, 4.0) );
	BOOST_CHECK_NE(    coord(8.0, 1.0, 5.0), coord(8.0, 1.0, 4.0) );
}

/// \brief Check multiplication
BOOST_AUTO_TEST_CASE(multiplication) {
	BOOST_CHECK_EQUAL( double_coord1, 2.0 * coord1  );
	BOOST_CHECK_EQUAL( double_coord1, coord1 * 2.0  );
	BOOST_CHECK_EQUAL( double_coord1, coord1 / 0.5  );
	BOOST_CHECK_NE(    double_coord1, 4.0 * coord1  );
	BOOST_CHECK_NE(    double_coord1, coord1 * 4.0  );
	BOOST_CHECK_NE(    double_coord1, coord1 / 0.25 );
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
	BOOST_CHECK_LT(                   dot_product(coord2,            coord_1_perp_to_2                ), ACCURACY_PERCENTAGE()/100.0 );
	BOOST_CHECK_LT(                   dot_product(coord_1_para_to_2, coord_1_perp_to_2                ), ACCURACY_PERCENTAGE()/100.0 );
	BOOST_CHECK_CLOSE(length(coord2), dot_product(coord2,            normalise_copy(coord_1_para_to_2)), ACCURACY_PERCENTAGE()       );
}

/// \brief Check dot product of orthogonal vectors
BOOST_AUTO_TEST_CASE(dot_products) {
	BOOST_CHECK_CLOSE(dot_product(coord(1.41421356237, 1.41421356237, 0.0), coord(2.0, 0.0, 0.0)), 2.82842712475, ACCURACY_PERCENTAGE());
	BOOST_CHECK_EQUAL(dot_product(coord(1.0, 0.0, 0.0), coord(1.0, 0.0, 0.0)), 1.0);
	BOOST_CHECK_EQUAL(dot_product(coord(0.0, 1.0, 0.0), coord(0.0, 1.0, 0.0)), 1.0);
	BOOST_CHECK_EQUAL(dot_product(coord(0.0, 0.0, 1.0), coord(0.0, 0.0, 1.0)), 1.0);
	BOOST_CHECK_LT(dot_product(coord1, coord2), ACCURACY_PERCENTAGE()/100.0);
	BOOST_CHECK_LT(dot_product(coord2, coord3), ACCURACY_PERCENTAGE()/100.0);
	BOOST_CHECK_LT(dot_product(coord3, coord1), ACCURACY_PERCENTAGE()/100.0);
}

/// \brief Check cross_product
BOOST_AUTO_TEST_CASE(cross_products) {
	BOOST_CHECK_CLOSE(cross_product(coord1, coord2).get_x(), coord3.get_x(), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(cross_product(coord1, coord2).get_y(), coord3.get_y(), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(cross_product(coord1, coord2).get_z(), coord3.get_z(), ACCURACY_PERCENTAGE());

	BOOST_CHECK_CLOSE(cross_product(coord2, coord3).get_x(), coord1.get_x(), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(cross_product(coord2, coord3).get_y(), coord1.get_y(), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(cross_product(coord2, coord3).get_z(), coord1.get_z(), ACCURACY_PERCENTAGE());

	BOOST_CHECK_CLOSE(cross_product(coord3, coord1).get_x(), coord2.get_x(), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(cross_product(coord3, coord1).get_y(), coord2.get_y(), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(cross_product(coord3, coord1).get_z(), coord2.get_z(), ACCURACY_PERCENTAGE());

	BOOST_CHECK_CLOSE(length(cross_product(coord1, coord2)), length(coord3), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(length(cross_product(coord2, coord3)), length(coord1), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(length(cross_product(coord3, coord1)), length(coord2), ACCURACY_PERCENTAGE());

	BOOST_CHECK_CLOSE(length(cross_product(2.0 * coord1, 3.0 * coord2)), length(6.0 * coord3), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(length(cross_product(2.0 * coord2, 3.0 * coord3)), length(6.0 * coord1), ACCURACY_PERCENTAGE());
	BOOST_CHECK_CLOSE(length(cross_product(2.0 * coord3, 3.0 * coord1)), length(6.0 * coord2), ACCURACY_PERCENTAGE());
}
/// \brief Check angle_between_two_vectors() works as expected
BOOST_AUTO_TEST_CASE(angle_between_two_vectors_works) {
	using coord_coord_pair = pair<coord, coord>;
	using coord_coord_vec = vector<coord_coord_pair>;
	const coord_coord_vec coord_pairs = {
		{ coord(  5,  0,  0 ), coord( -5,  0,  0 ) },
		{ coord(  0,  5,  0 ), coord(  0, -5,  0 ) },
		{ coord(  0,  0,  5 ), coord(  0,  0, -5 ) }
	};
	for (const coord_coord_pair &coord_pair_1 : coord_pairs) {
		BOOST_CHECK_EQUAL( one_revolution<double>() / 2.0, angle_between_two_vectors( coord_pair_1.first,  coord_pair_1.second ) );
		BOOST_CHECK_EQUAL( one_revolution<double>() / 2.0, angle_between_two_vectors( coord_pair_1.second, coord_pair_1.first  ) );

		for (const coord_coord_pair &coord_pair_2 : coord_pairs) {
			if ( coord_pair_1.first != coord_pair_2.first ) {
				BOOST_CHECK_EQUAL( one_revolution<double>() / 4.0, angle_between_two_vectors( coord_pair_1.first,  coord_pair_2.first  ) );
				BOOST_CHECK_EQUAL( one_revolution<double>() / 4.0, angle_between_two_vectors( coord_pair_1.first,  coord_pair_2.second ) );
				BOOST_CHECK_EQUAL( one_revolution<double>() / 4.0, angle_between_two_vectors( coord_pair_1.second, coord_pair_2.first  ) );
				BOOST_CHECK_EQUAL( one_revolution<double>() / 4.0, angle_between_two_vectors( coord_pair_1.second, coord_pair_2.second ) );
			}
		}
	}
}

/// \brief Check dihedral_angle_between_four_points() works as expected
BOOST_AUTO_TEST_CASE(dihedral_angle_between_four_points_works) {
	BOOST_CHECK_CLOSE(
		angle_in_degrees(       - one_revolution<double>() / 8.0 ),
		angle_in_degrees( dihedral_angle_between_four_points(coord(0, -1,  1), coord(0, 0, 0), coord(1, 0, 0), coord(1, 0, 1) ) ),
		ACCURACY_PERCENTAGE()
	);
	BOOST_CHECK_CLOSE(
		angle_in_degrees(         one_revolution<double>() / 8.0 ),
		angle_in_degrees( dihedral_angle_between_four_points(coord(0,  1,  1), coord(0, 0, 0), coord(1, 0, 0), coord(1, 0, 1) ) ),
		ACCURACY_PERCENTAGE()
	);
	BOOST_CHECK_CLOSE(
		angle_in_degrees( - 3.0 * one_revolution<double>() / 8.0 ),
		angle_in_degrees( dihedral_angle_between_four_points(coord(0, -1, -1), coord(0, 0, 0), coord(1, 0, 0), coord(1, 0, 1) ) ),
		ACCURACY_PERCENTAGE()
	);
	BOOST_CHECK_CLOSE(
		angle_in_degrees(   3.0 * one_revolution<double>() / 8.0 ),
		angle_in_degrees( dihedral_angle_between_four_points(coord(0,  1, -1), coord(0, 0, 0), coord(1, 0, 0), coord(1, 0, 1) ) ),
		ACCURACY_PERCENTAGE()
	);
}

BOOST_AUTO_TEST_SUITE_END()
