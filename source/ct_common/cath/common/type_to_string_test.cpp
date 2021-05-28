/// \file
/// \brief The type_to_string test suite

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

#include "type_to_string.hpp"

#include <boost/test/unit_test.hpp>

using namespace ::cath::common;

using ::std::literals::string_literals::operator""s;
using ::std::tuple;
using ::std::vector;

BOOST_AUTO_TEST_SUITE(type_to_string_test_suite)

BOOST_AUTO_TEST_CASE(int_works) {
	BOOST_CHECK_EQUAL( type_to_string<int>(), "int"s );
}

BOOST_AUTO_TEST_CASE(int_vec_works) {
	BOOST_CHECK_EQUAL( type_to_string<vector<int>>(), "std::vector<int, std::allocator<int>>"s );
}

BOOST_AUTO_TEST_CASE(int_doub_float_tuple_works) {
	using int_doub_float_tuple = tuple<int, double, float>;
	BOOST_CHECK_EQUAL( type_to_string<int_doub_float_tuple>(), "std::tuple<int, double, float>"s );
}


BOOST_AUTO_TEST_SUITE_END()
