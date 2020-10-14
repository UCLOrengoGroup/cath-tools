/// \file
/// \brief The transform_build test suite

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

#include "is_tuple.hpp"

#include <boost/test/unit_test.hpp>

using namespace ::cath::common;

using ::std::tuple;
using ::std::vector;

BOOST_AUTO_TEST_SUITE(is_tuple_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert( ! is_tuple            <               int     >::value, "              int     should not be identified as a tuple" );
	static_assert( ! is_tuple            < const         int     >::value, "const         int     should not be identified as a tuple" );
	static_assert( ! is_tuple            <       vector< int >   >::value, "      vector< int >   should not be identified as a tuple" );
	static_assert( ! is_tuple            < const vector< int >   >::value, "const vector< int >   should not be identified as a tuple" );

	static_assert(   is_tuple            <       tuple < int >   >::value, "      tuple < int >   should     be identified as a tuple" );
	static_assert( ! is_tuple            < const tuple < int >   >::value, "const tuple < int >   should not be identified as a tuple" );
	static_assert( ! is_tuple            <       tuple < int > & >::value, "      tuple < int > & should not be identified as a tuple" );
	static_assert( ! is_tuple            < const tuple < int > & >::value, "const tuple < int > & should not be identified as a tuple" );

	static_assert( ! is_tuple_after_decay<               int     >::value, "              int     should not be identified as a tuple" );
	static_assert( ! is_tuple_after_decay< const         int     >::value, "const         int     should not be identified as a tuple" );
	static_assert( ! is_tuple_after_decay<       vector< int >   >::value, "      vector< int >   should not be identified as a tuple" );
	static_assert( ! is_tuple_after_decay< const vector< int >   >::value, "const vector< int >   should not be identified as a tuple" );

	static_assert(   is_tuple_after_decay<       tuple < int >   >::value, "      tuple < int >   should     be identified as a tuple" );
	static_assert(   is_tuple_after_decay< const tuple < int >   >::value, "const tuple < int >   should     be identified as a tuple" );
	static_assert(   is_tuple_after_decay<       tuple < int > & >::value, "      tuple < int > & should     be identified as a tuple" );
	static_assert(   is_tuple_after_decay< const tuple < int > & >::value, "const tuple < int > & should     be identified as a tuple" );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
