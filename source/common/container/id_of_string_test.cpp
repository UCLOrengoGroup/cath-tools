/// \file
/// \brief The id_of_string test suite

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

#include <boost/range/irange.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/container/id_of_string.hpp"
#include "common/size_t_literal.hpp"

namespace cath { namespace test { } }

using namespace cath::common;
using namespace cath::test;

using boost::irange;

namespace cath {
	namespace test {

		/// \brief The id_of_string_test_suite_fixture to assist in testing id_of_string
		struct id_of_string_test_suite_fixture {
		protected:
			~id_of_string_test_suite_fixture() noexcept = default;

		public:
			const std::string example_string_zero  = "this is zero"  ;
			const std::string example_string_one   = "this is one"   ;
			const std::string example_string_two   = "this is two"   ;
			const std::string example_string_three = "this is three" ;
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(id_of_string_test_suite, id_of_string_test_suite_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	id_of_string the_ios;

	BOOST_CHECK      (   the_ios.empty()   );
	BOOST_CHECK_EQUAL(   the_ios.size(), 0 );

	BOOST_CHECK_EQUAL(   the_ios.emplace( example_string_zero ), 0 );

	BOOST_CHECK      ( ! the_ios.empty()   );
	BOOST_CHECK_EQUAL(   the_ios.size(), 1 );
	BOOST_CHECK_EQUAL(   the_ios.get( example_string_zero ), 0 );

	BOOST_CHECK_EQUAL(   the_ios.emplace_if_not_present( example_string_zero ), 0 );

	BOOST_CHECK      ( ! the_ios.empty()   );
	BOOST_CHECK_EQUAL(   the_ios.size(), 1 );
	BOOST_CHECK_EQUAL(   the_ios.get( example_string_zero ), 0 );

	BOOST_CHECK_EQUAL(   the_ios.emplace_if_not_present( example_string_one ), 1 );

	BOOST_CHECK      ( ! the_ios.empty()   );
	BOOST_CHECK_EQUAL(   the_ios.size(), 2 );
	BOOST_CHECK_EQUAL(   the_ios.get( example_string_zero ), 0 );
	BOOST_CHECK_EQUAL(   the_ios.get( example_string_one  ), 1 );

	BOOST_CHECK_EQUAL(   the_ios.emplace( example_string_two   ), 2 );
	BOOST_CHECK_EQUAL(   the_ios.emplace( example_string_three ), 3 );
}

BOOST_AUTO_TEST_CASE(fast) {
	using std::to_string;
	constexpr size_t MAX = 1000_z;

	id_of_string the_ios;
	the_ios.reserve( MAX );
	for (const size_t &ctr : irange( 0_z, MAX ) ) {
		the_ios.emplace( to_string( ctr ) );
	}
	BOOST_CHECK_EQUAL( the_ios.get(   "7" ),   7 );
	BOOST_CHECK_EQUAL( the_ios.get( "997" ), 997 );
}

BOOST_AUTO_TEST_SUITE_END()

