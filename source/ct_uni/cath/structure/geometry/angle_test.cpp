/// \file
/// \brief The angle test suite

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

#include <boost/math/constants/constants.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/test/unit_test.hpp>

//#include <boost/units/systems/si/plane_angle.hpp>

#include "cath/common/difference.hpp"
#include "cath/structure/geometry/angle.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "cath/test/global_test_constants.hpp"
#include "cath/test/test_tools.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::common::test;
using namespace ::cath::geom;

using ::boost::math::constants::pi;
using ::std::array;

/// \brief Test suite for the angle
BOOST_FIXTURE_TEST_SUITE(angle_test_suite, global_test_constants)

using angle_value_types = boost::mpl::vector<double, float>;

/// \brief Check that the ctor does not throw
BOOST_AUTO_TEST_CASE_TEMPLATE(ctor_does_not_throw, T, angle_value_types) {
	[[maybe_unused]] constexpr auto a = angle<T>( 1.0 );
	BOOST_TEST (true );
}

/// \brief Check that the factories do not throw
BOOST_AUTO_TEST_CASE_TEMPLATE(factories_do_not_throw, T, angle_value_types) {
	[[maybe_unused]] constexpr auto a1 = make_angle_from_degrees    <T>( 1.0 );
	[[maybe_unused]] constexpr auto a2 = make_angle_from_radians    <T>( 1.0 );
	[[maybe_unused]] constexpr auto a3 = make_angle_from_revolutions<T>( 1.0 );
	BOOST_TEST (true );
}

/// \brief Check that equality/inequality works
BOOST_AUTO_TEST_CASE_TEMPLATE(equality_works, T, angle_value_types) {
	constexpr array angles = {
		make_angle_from_degrees    <T>( 1.0 ),
		make_angle_from_radians    <T>( 1.0 ),
		make_angle_from_revolutions<T>( 1.0 )
	};
	check_equality_operators_on_diff_vals_range( angles );
}

/// \brief Check that sensible things are equal
BOOST_AUTO_TEST_CASE_TEMPLATE(angle_values, T, angle_value_types) {
	static_assert( make_angle_from_degrees    <T>(                 0.0 ) == ZERO_ANGLE<T>     );
	static_assert( make_angle_from_degrees    <T>(               360.0 ) == ONE_REVOLUTION<T> );

	static_assert( make_angle_from_degrees    <T>(              -360.0 ) == make_angle_from_radians    <T>( pi<double>() * -2.0 ) );
	static_assert( make_angle_from_radians    <T>( pi<double>() * -2.0 ) == make_angle_from_revolutions<T>(                -1.0 ) );
	static_assert( make_angle_from_revolutions<T>(                -1.0 ) == make_angle_from_degrees    <T>(              -360.0 ) );

	static_assert( make_angle_from_degrees    <T>(              -180.0 ) == make_angle_from_radians    <T>( pi<double>() * -1.0 ) );
	static_assert( make_angle_from_radians    <T>( pi<double>() * -1.0 ) == make_angle_from_revolutions<T>(                -0.5 ) );
	static_assert( make_angle_from_revolutions<T>(                -0.5 ) == make_angle_from_degrees    <T>(              -180.0 ) );

	static_assert( make_angle_from_degrees    <T>(                 0.0 ) == make_angle_from_radians    <T>(                 0.0 ) );
	static_assert( make_angle_from_radians    <T>(                 0.0 ) == make_angle_from_revolutions<T>(                 0.0 ) );
	static_assert( make_angle_from_revolutions<T>(                 0.0 ) == make_angle_from_degrees    <T>(                 0.0 ) );

	static_assert( make_angle_from_degrees    <T>(               180.0 ) == make_angle_from_radians    <T>(  pi<double>() * 1.0 ) );
	static_assert( make_angle_from_radians    <T>(  pi<double>() * 1.0 ) == make_angle_from_revolutions<T>(                 0.5 ) );
	static_assert( make_angle_from_revolutions<T>(                 0.5 ) == make_angle_from_degrees    <T>(               180.0 ) );

	static_assert( make_angle_from_degrees    <T>(               360.0 ) == make_angle_from_radians    <T>(  pi<double>() * 2.0 ) );
	static_assert( make_angle_from_radians    <T>(  pi<double>() * 2.0 ) == make_angle_from_revolutions<T>(                 1.0 ) );
	static_assert( make_angle_from_revolutions<T>(                 1.0 ) == make_angle_from_degrees    <T>(               360.0 ) );

	BOOST_TEST (true );
}

/// \brief Check that sensible things are equal
BOOST_AUTO_TEST_CASE_TEMPLATE(shift_works, T, angle_value_types) {
	static_assert( shift_copy(     ZERO_ANGLE<T>                                                ) ==  ZERO_ANGLE<T>     );
	static_assert( shift_copy(     ZERO_ANGLE<T>, ZERO_ANGLE<T>                                 ) ==  ZERO_ANGLE<T>     );
	static_assert( shift_copy(     ZERO_ANGLE<T>, ZERO_ANGLE<T>, angle_endpoint_loc::USE_LOWER  ) ==  ZERO_ANGLE<T>     );
	static_assert( shift_copy(     ZERO_ANGLE<T>, ZERO_ANGLE<T>, angle_endpoint_loc::USE_UPPER  ) ==  ONE_REVOLUTION<T> );
	static_assert( shift_copy(     ZERO_ANGLE<T>, ZERO_ANGLE<T>, angle_endpoint_loc::USE_EITHER ) ==  ZERO_ANGLE<T>     );

	static_assert( shift_copy( ONE_REVOLUTION<T>                                                ) ==  ZERO_ANGLE<T>     );
	static_assert( shift_copy( ONE_REVOLUTION<T>, ZERO_ANGLE<T>                                 ) ==  ZERO_ANGLE<T>     );
	static_assert( shift_copy( ONE_REVOLUTION<T>, ZERO_ANGLE<T>, angle_endpoint_loc::USE_LOWER  ) ==  ZERO_ANGLE<T>     );
	static_assert( shift_copy( ONE_REVOLUTION<T>, ZERO_ANGLE<T>, angle_endpoint_loc::USE_UPPER  ) ==  ONE_REVOLUTION<T> );

//	cerr << "*********** About to enter troublesome case ***********" << endl;
	static_assert( shift_copy( ONE_REVOLUTION<T>, ZERO_ANGLE<T>, angle_endpoint_loc::USE_EITHER ) ==  ONE_REVOLUTION<T> );
//	cerr << "*********** Finished troublesome case ***********" << endl;

	static_assert( shift_copy( make_angle_from_revolutions<double>( -0.5 )                      ) ==  make_angle_from_revolutions<double>( 0.5) );
	static_assert( shift_copy( make_angle_from_revolutions<double>(  1.5 )                      ) ==  make_angle_from_revolutions<double>( 0.5) );

	BOOST_TEST (true );
}

/// \brief Check that sensible things are equal
BOOST_AUTO_TEST_CASE_TEMPLATE(difference_and_wrapped_difference, T, angle_value_types) {
	constexpr auto one_quarter_turn    = make_angle_from_revolutions<T>( 0.25  );
	constexpr auto three_quarters_turn = make_angle_from_revolutions<T>( 0.75  );

	constexpr auto one_eighth_turn     = make_angle_from_revolutions<T>( 0.875 );
	constexpr auto seven_eighths_turn  = make_angle_from_revolutions<T>( 0.125 );

	static_assert( difference( seven_eighths_turn,    one_eighth_turn ) == three_quarters_turn );
	static_assert( difference(    one_eighth_turn, seven_eighths_turn ) == three_quarters_turn );
	BOOST_CHECK_CLOSE( angle_in_degrees( wrapped_difference( seven_eighths_turn,    one_eighth_turn ) ), angle_in_degrees(    one_quarter_turn ), LOOSER_ACCURACY_PERCENTAGE );
	BOOST_CHECK_CLOSE( angle_in_degrees( wrapped_difference(    one_eighth_turn, seven_eighths_turn ) ), angle_in_degrees(    one_quarter_turn ), LOOSER_ACCURACY_PERCENTAGE );
}

BOOST_AUTO_TEST_SUITE_END()
