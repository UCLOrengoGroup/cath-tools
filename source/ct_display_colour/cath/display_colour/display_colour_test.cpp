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

#include <array>
#include <tuple>

#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/boost_addenda/range/to_vector.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/display_colour/display_colour.hpp"
#include "cath/display_colour/display_colour_type_aliases.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "cath/test/global_test_constants.hpp"
#include "cath/test/test_tools.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::common::test;

using ::boost::lexical_cast;
using ::std::apply;
using ::std::array;
using ::std::string;

namespace {

	namespace detail {

		/// Make a copy of the specified array with Extra extra slots at the end
		template <size_t Extra, typename T, size_t N>
		constexpr array<T, N + Extra> with_n_more_at_end( const ::std::array<T, N> &prm_array ) {
			return apply( []( const auto &...x ) { return array<T, N + Extra>{ x... }; }, prm_array );
		}

		constexpr array<double, global_test_constants::INVALID_DOUBLES.size() + 2> make_invalid_components() {
			auto invalid_components = with_n_more_at_end<2>( global_test_constants::INVALID_DOUBLES );

			// -0.1 should be invalid because colour components must be >= 0.0
			invalid_components.at( global_test_constants::INVALID_DOUBLES.size() ) = -0.1;
			//  1.1 should be invalid because colour components must be <= 1.0
			invalid_components.at( global_test_constants::INVALID_DOUBLES.size() + 1 ) = 1.1;
			return invalid_components;
		}

	} // namespace detail

	/// \brief The display_colour_test_suite_fixture to assist in testing display_colour_test
	struct display_colour_test_suite_fixture : protected global_test_constants {
	  protected:
		~display_colour_test_suite_fixture() noexcept = default;

	  public:
		/// \brief Check that the r, g and b component values of a viewer colour are as expected
		void check_r_g_b_of_display_colour( const display_colour &prm_display_colour, ///< The display_colour to check
		                                    const double &        prm_expected_r, ///< The expected r component value
		                                    const double &        prm_expected_g, ///< The expected g component value
		                                    const double &        prm_expected_b  ///< The expected b component value
		                                    ) {
			BOOST_CHECK_EQUAL( prm_expected_r, prm_display_colour.get_r() );
			BOOST_CHECK_EQUAL( prm_expected_g, prm_display_colour.get_g() );
			BOOST_CHECK_EQUAL( prm_expected_b, prm_display_colour.get_b() );
		}

		static constexpr double EXPECTED_R   = 0.684;
		static constexpr double EXPECTED_G   = 0.334;
		static constexpr double EXPECTED_B   = 0.638;
		static constexpr array  EXPECTED_RGB = { EXPECTED_R, EXPECTED_G, EXPECTED_B };

		static constexpr auto INVALID_COMPONENT_VALUES = detail::make_invalid_components();

		const string EXPECTED_RGB_STRING = { lexical_cast<string>( EXPECTED_R ) + "," + lexical_cast<string>( EXPECTED_G )
			                                 + "," + lexical_cast<string>( EXPECTED_B ) };
	};

} // namespace

BOOST_FIXTURE_TEST_SUITE(display_colour_test_suite, display_colour_test_suite_fixture)

/// \brief Check that the ctor throws an invalid_argument_exception if given any invalid value
BOOST_AUTO_TEST_CASE(ctor_rejects_invalid_values) {
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
	const auto standard_colours = { RED,
	                                YELLOW,
	                                GREEN,
	                                CYAN,
	                                BLUE,
	                                MAGENTA };
	check_equality_operators_on_diff_vals_range( standard_colours );
}

BOOST_AUTO_TEST_SUITE_END()

