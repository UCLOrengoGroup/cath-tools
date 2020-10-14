/// \file
/// \brief The booled_to_string test suite

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#include <boost/mpl/vector.hpp>
#include <boost/range/irange.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/string/booled_to_string.hpp"

namespace cath { namespace test {} }

using namespace cath::common;
using namespace cath::test;

using boost::irange;
using std::to_string;

namespace cath {
	namespace test {
		/// \brief A bunch of type for which booled_to_string should handle
		using booled_to_string_int_types = boost::mpl::vector<char,
		                                                      short int,
		                                                      short unsigned int,
		                                                      int,
		                                                      unsigned int,
		                                                      long int,
		                                                      long unsigned int>;
	}
}  // namespace cath

BOOST_AUTO_TEST_SUITE(booled_to_string_test_suite)

/// \brief Check that booled_to_string() correctly returns "true" for bool true
BOOST_AUTO_TEST_CASE(works_for_bool_true) {
	BOOST_CHECK_EQUAL( booled_to_string( true ), "true" );
}

/// \brief Check that booled_to_string() correctly returns "false" for bool false
BOOST_AUTO_TEST_CASE(works_for_bool_false) {
	BOOST_CHECK_EQUAL( booled_to_string( false  ), "false"  );
}

/// \brief Check that booled_to_string() correctly returns "0" for int_type 0
BOOST_AUTO_TEST_CASE_TEMPLATE(works_for_int_type_0, int_type, booled_to_string_int_types) {
	BOOST_CHECK_EQUAL( booled_to_string( static_cast<int_type>( 0 ) ), "0" );
}

/// \brief Check that booled_to_string() correctly returns "1" for int_type 1
BOOST_AUTO_TEST_CASE_TEMPLATE(works_for_int_type_1, int_type, booled_to_string_int_types) {
	BOOST_CHECK_EQUAL( booled_to_string( static_cast<int_type>( 1 ) ), "1" );
}

/// \brief Check that booled_to_string() correctly returns the same as to_string() for int_types in 0 to 10 inclusive
BOOST_AUTO_TEST_CASE_TEMPLATE(works_for_int_types_0_to_10, int_type, booled_to_string_int_types) {
	for (const int_type &x : irange( static_cast<int_type>( 0 ), static_cast<int_type>( 11 ) ) ) {
		BOOST_CHECK_EQUAL(
			booled_to_string( static_cast<int_type>( x ) ),
			to_string       ( static_cast<int_type>( x ) )
		);
	}

}

BOOST_AUTO_TEST_SUITE_END()
