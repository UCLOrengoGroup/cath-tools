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

#include <boost/math/constants/constants.hpp>
#include <boost/test/auto_unit_test.hpp>

//#include <boost/units/systems/si/plane_angle.hpp>

#include "common/boost_addenda/test/boost_check_no_throw_diag.hpp"
#include "common/difference.hpp"
#include "common/test_tools.hpp"
#include "structure/geometry/angle.hpp"
#include "test/global_test_constants.hpp"

using namespace boost::math::constants;
using namespace cath;
using namespace cath::common;
using namespace cath::common::test;
using namespace cath::geom;

namespace cath {
	namespace test {

		/// \brief The angle_test_suite_fixture to assist in testing test
		struct angle_test_suite_fixture : protected global_test_constants {
		protected:
			~angle_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

/// \brief Test suite for the angle
BOOST_FIXTURE_TEST_SUITE(angle_test_suite, cath::test::angle_test_suite_fixture)

using angle_value_types = boost::mpl::vector<double, float>;

/// \brief Check that the ctor does not throw
BOOST_AUTO_TEST_CASE_TEMPLATE(ctor_does_not_throw, T, angle_value_types) {
//	BOOST_CHECK_NO_THROW_DIAG( angle the_angle( 1.0 * boost::units::si::radian ) );
	BOOST_CHECK_NO_THROW_DIAG( angle<T> the_angle( 1.0 ) );
}

/// \brief Check that the factories do not throw
BOOST_AUTO_TEST_CASE_TEMPLATE(factories_do_not_throw, T, angle_value_types) {
	BOOST_CHECK_NO_THROW_DIAG( make_angle_from_degrees    <T>( 1.0 ) );
	BOOST_CHECK_NO_THROW_DIAG( make_angle_from_radians    <T>( 1.0 ) );
	BOOST_CHECK_NO_THROW_DIAG( make_angle_from_revolutions<T>( 1.0 ) );
}

/// \brief Check that equality/inequality works
BOOST_AUTO_TEST_CASE_TEMPLATE(equality_works, T, angle_value_types) {
	const auto angles = {
		make_angle_from_degrees    <T>( 1.0 ),
		make_angle_from_radians    <T>( 1.0 ),
		make_angle_from_revolutions<T>( 1.0 )
	};
	check_equality_operators_on_diff_vals_range( angles );
}

/// \brief Check that sensible things are equal
BOOST_AUTO_TEST_CASE_TEMPLATE(angle_values, T, angle_value_types) {
	BOOST_CHECK_EQUAL( make_angle_from_degrees    <T>(                 0.0 ), zero_angle<T>()     );
	BOOST_CHECK_EQUAL( make_angle_from_degrees    <T>(               360.0 ), one_revolution<T>() );

	BOOST_CHECK_EQUAL( make_angle_from_degrees    <T>(              -360.0 ), make_angle_from_radians    <T>( pi<double>() * -2.0 ) );
	BOOST_CHECK_EQUAL( make_angle_from_radians    <T>( pi<double>() * -2.0 ), make_angle_from_revolutions<T>(                -1.0 ) );
	BOOST_CHECK_EQUAL( make_angle_from_revolutions<T>(                -1.0 ), make_angle_from_degrees    <T>(              -360.0 ) );

	BOOST_CHECK_EQUAL( make_angle_from_degrees    <T>(              -180.0 ), make_angle_from_radians    <T>( pi<double>() * -1.0 ) );
	BOOST_CHECK_EQUAL( make_angle_from_radians    <T>( pi<double>() * -1.0 ), make_angle_from_revolutions<T>(                -0.5 ) );
	BOOST_CHECK_EQUAL( make_angle_from_revolutions<T>(                -0.5 ), make_angle_from_degrees    <T>(              -180.0 ) );

	BOOST_CHECK_EQUAL( make_angle_from_degrees    <T>(                 0.0 ), make_angle_from_radians    <T>(                 0.0 ) );
	BOOST_CHECK_EQUAL( make_angle_from_radians    <T>(                 0.0 ), make_angle_from_revolutions<T>(                 0.0 ) );
	BOOST_CHECK_EQUAL( make_angle_from_revolutions<T>(                 0.0 ), make_angle_from_degrees    <T>(                 0.0 ) );

	BOOST_CHECK_EQUAL( make_angle_from_degrees    <T>(               180.0 ), make_angle_from_radians    <T>(  pi<double>() * 1.0 ) );
	BOOST_CHECK_EQUAL( make_angle_from_radians    <T>(  pi<double>() * 1.0 ), make_angle_from_revolutions<T>(                 0.5 ) );
	BOOST_CHECK_EQUAL( make_angle_from_revolutions<T>(                 0.5 ), make_angle_from_degrees    <T>(               180.0 ) );

	BOOST_CHECK_EQUAL( make_angle_from_degrees    <T>(               360.0 ), make_angle_from_radians    <T>(  pi<double>() * 2.0 ) );
	BOOST_CHECK_EQUAL( make_angle_from_radians    <T>(  pi<double>() * 2.0 ), make_angle_from_revolutions<T>(                 1.0 ) );
	BOOST_CHECK_EQUAL( make_angle_from_revolutions<T>(                 1.0 ), make_angle_from_degrees    <T>(               360.0 ) );
}

/// \brief Check that sensible things are equal
BOOST_AUTO_TEST_CASE_TEMPLATE(shift_works, T, angle_value_types) {
	BOOST_CHECK_EQUAL( shift_copy(     zero_angle<T>()                                                  ), zero_angle<T>()     );
	BOOST_CHECK_EQUAL( shift_copy(     zero_angle<T>(), zero_angle<T>()                                 ), zero_angle<T>()     );
	BOOST_CHECK_EQUAL( shift_copy(     zero_angle<T>(), zero_angle<T>(), angle_endpoint_loc::USE_LOWER  ), zero_angle<T>()     );
	BOOST_CHECK_EQUAL( shift_copy(     zero_angle<T>(), zero_angle<T>(), angle_endpoint_loc::USE_UPPER  ), one_revolution<T>() );
	BOOST_CHECK_EQUAL( shift_copy(     zero_angle<T>(), zero_angle<T>(), angle_endpoint_loc::USE_EITHER ), zero_angle<T>()     );

	BOOST_CHECK_EQUAL( shift_copy( one_revolution<T>()                                                  ), zero_angle<T>()     );
	BOOST_CHECK_EQUAL( shift_copy( one_revolution<T>(), zero_angle<T>()                                 ), zero_angle<T>()     );
	BOOST_CHECK_EQUAL( shift_copy( one_revolution<T>(), zero_angle<T>(), angle_endpoint_loc::USE_LOWER  ), zero_angle<T>()     );
	BOOST_CHECK_EQUAL( shift_copy( one_revolution<T>(), zero_angle<T>(), angle_endpoint_loc::USE_UPPER  ), one_revolution<T>() );

//	cerr << "*********** About to enter troublesome case ***********" << endl;
	BOOST_CHECK_EQUAL( shift_copy( one_revolution<T>(), zero_angle<T>(), angle_endpoint_loc::USE_EITHER ), one_revolution<T>() );
//	cerr << "*********** Finished troublesome case ***********" << endl;

	BOOST_CHECK_EQUAL( shift_copy( make_angle_from_revolutions<double>( -0.5 )                 ), make_angle_from_revolutions<double>( 0.5) );
	BOOST_CHECK_EQUAL( shift_copy( make_angle_from_revolutions<double>(  1.5 )                 ), make_angle_from_revolutions<double>( 0.5) );
}

/// \brief Check that sensible things are equal
BOOST_AUTO_TEST_CASE_TEMPLATE(difference_and_wrapped_difference, T, angle_value_types) {
	const auto one_quarter_turn    = make_angle_from_revolutions<T>( 0.25  );
	const auto three_quarters_turn = make_angle_from_revolutions<T>( 0.75  );

	const auto one_eighth_turn     = make_angle_from_revolutions<T>( 0.875 );
	const auto seven_eighths_turn  = make_angle_from_revolutions<T>( 0.125 );

	BOOST_CHECK_EQUAL(                           difference( seven_eighths_turn,    one_eighth_turn ),                     three_quarters_turn                          );
	BOOST_CHECK_EQUAL(                           difference(    one_eighth_turn, seven_eighths_turn ),                     three_quarters_turn                          );
	BOOST_CHECK_CLOSE( angle_in_degrees( wrapped_difference( seven_eighths_turn,    one_eighth_turn ) ), angle_in_degrees(    one_quarter_turn ), LOOSER_ACCURACY_PERCENTAGE() );
	BOOST_CHECK_CLOSE( angle_in_degrees( wrapped_difference(    one_eighth_turn, seven_eighths_turn ) ), angle_in_degrees(    one_quarter_turn ), LOOSER_ACCURACY_PERCENTAGE() );
}

BOOST_AUTO_TEST_SUITE_END()
