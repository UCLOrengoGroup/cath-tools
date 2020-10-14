/// \file
/// \brief The are_same test suite

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

#include "are_same.hpp"

#include <boost/test/unit_test.hpp>

using namespace ::cath::common;

using ::std::equal_to;
using ::std::vector;

BOOST_AUTO_TEST_SUITE(are_same_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	vector<int> a = {  2,  2,  2, -2,  2 };
	vector<int> b = {  2,  2,  2,  2,  2 };

	BOOST_TEST( ! are_same( a, equal_to<>{} ) );
	BOOST_TEST( ! are_same( a               ) );

	BOOST_TEST(   are_same( b, equal_to<>{} ) );
	BOOST_TEST(   are_same( b               ) );

	BOOST_TEST(   are_same( a, equal_to<>{}, [] (const int &x) { return x * x; } ) );
	BOOST_TEST( ! are_same( a, equal_to<>{}, [] (const int &x) { return x + x; } ) );
}

BOOST_AUTO_TEST_SUITE_END()
