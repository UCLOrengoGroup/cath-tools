/// \file
/// \brief The display_colour_test test suite

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

#include <boost/lexical_cast.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/boost_check_no_throw_diag.h"
#include "common/test_tools.h"
#include "exception/invalid_argument_exception.h"
#include "test/global_test_constants.h"
#include "display/display_type_aliases.h"
#include "display/display_colour/display_colour.h"

using namespace cath;
using namespace cath::common;
using namespace cath::common::test;
using namespace std;

using boost::lexical_cast;

namespace cath {
	namespace test {

		/// \brief The display_colour_test_suite_fixture to assist in testing display_colour_test
		struct display_colour_test_suite_fixture : public global_test_constants {
		private:
			doub_vec get_invalid_component_values() const {
				// The standard invalid doubles
				doub_vec invalid_doubles = INVALID_DOUBLES();
				// -0.1 should be invalid because colour components must be >= 0.0
				invalid_doubles.push_back( -0.1 );
				//  1.1 should be invalid because colour components must be <= 1.0
				return invalid_doubles;
			}

		protected:
			~display_colour_test_suite_fixture() noexcept = default;

		public:
			/// \brief Check that the r, g and b component values of a viewer colour are as expected
			void check_r_g_b_of_display_colour(const display_colour &arg_display_colour, ///< The display_colour to check
			                                  const double        &arg_expected_r,    ///< The expected r component value
			                                  const double        &arg_expected_g,    ///< The expected g component value
			                                  const double        &arg_expected_b     ///< The expected b component value
			                                  ) {
				BOOST_CHECK_EQUAL( arg_expected_r, arg_display_colour.get_r() );
				BOOST_CHECK_EQUAL( arg_expected_g, arg_display_colour.get_g() );
				BOOST_CHECK_EQUAL( arg_expected_b, arg_display_colour.get_b() );
			}

			const double   EXPECTED_R               = { 0.684 };
			const double   EXPECTED_G               = { 0.334 };
			const double   EXPECTED_B               = { 0.638 };
			const doub_vec EXPECTED_RGB             = { EXPECTED_R, EXPECTED_G, EXPECTED_B };
			const doub_vec INVALID_COMPONENT_VALUES = { get_invalid_component_values() };
			const string   EXPECTED_RGB_STRING      = {
				lexical_cast<string>( EXPECTED_R ) + "," +
				lexical_cast<string>( EXPECTED_G ) + "," +
				lexical_cast<string>( EXPECTED_B )
			};
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(display_colour_test_suite, cath::test::display_colour_test_suite_fixture)

/// \brief Check that the ctor throws an invalid_argument_exception if given any invalid value
BOOST_AUTO_TEST_CASE(ctor_rejects_invalid_values) {
	doub_vec invalid_doubles = INVALID_DOUBLES();
	invalid_doubles.push_back(-0.1);
	invalid_doubles.push_back( 0.1);
	for (const double &invalid_value : INVALID_COMPONENT_VALUES) {
		BOOST_CHECK_THROW(display_colour( invalid_value, EXPECTED_G,    EXPECTED_B    ), invalid_argument_exception);
		BOOST_CHECK_THROW(display_colour( EXPECTED_R,    invalid_value, EXPECTED_B    ), invalid_argument_exception);
		BOOST_CHECK_THROW(display_colour( EXPECTED_R,    EXPECTED_G,    invalid_value ), invalid_argument_exception);
	}
}

/// \brief Check that the ctor accepts valid values (at both extremes and in between)
BOOST_AUTO_TEST_CASE(ctor_accepts_valid_values) {
	BOOST_CHECK_NO_THROW_DIAG(display_colour(0.0, 0.5, 1.0));
	BOOST_CHECK_NO_THROW_DIAG(display_colour(1.0, 0.0, 0.5));
	BOOST_CHECK_NO_THROW_DIAG(display_colour(0.5, 1.0, 0.0));
}

/// \brief Check that the getters are const methods that return the correct values
BOOST_AUTO_TEST_CASE(getters_get_correct_values) {
	// Construct a const test_colour to check the getters are const methods
	const display_colour test_colour(
		EXPECTED_R,
		EXPECTED_G,
		EXPECTED_B
	);
	check_r_g_b_of_display_colour(
		test_colour,
		EXPECTED_R,
		EXPECTED_G,
		EXPECTED_B
	);
}

/// \brief Check that display_colour_from_components() works as expected from a const vector of three doubles
BOOST_AUTO_TEST_CASE(display_colour_from_components_works) {
	check_r_g_b_of_display_colour(
		display_colour_from_components(EXPECTED_RGB),
		EXPECTED_R,
		EXPECTED_G,
		EXPECTED_B
	);
}

/// \brief Check that display_colour_from_string() works as expected from a string of comma separated numbers
BOOST_AUTO_TEST_CASE(display_colour_from_string_works) {
	check_r_g_b_of_display_colour(
		display_colour_from_string(EXPECTED_RGB_STRING),
		EXPECTED_R,
		EXPECTED_G,
		EXPECTED_B
	);
}

/// \brief Check that comma_separated_string_of_display_colour() works by cheCking it reverses display_colour_from_string_works()
///        (which has its own independent test, above)
BOOST_AUTO_TEST_CASE(comma_separated_string_of_display_colour_works) {
	BOOST_CHECK_EQUAL(
		EXPECTED_RGB_STRING,
		comma_separated_string_of_display_colour(
			display_colour_from_string(
				EXPECTED_RGB_STRING
			)
		)
	);
}

/// \brief Check that the equality and inequality operators work as expected
BOOST_AUTO_TEST_CASE(equality_and_inequality_operators) {
	const auto standard_colours = { display_colour::RED,
	                                display_colour::YELLOW,
	                                display_colour::GREEN,
	                                display_colour::CYAN,
	                                display_colour::BLUE,
	                                display_colour::MAGENTA };
	check_equality_operators_on_diff_vals_range( standard_colours );
}

BOOST_AUTO_TEST_SUITE_END()

