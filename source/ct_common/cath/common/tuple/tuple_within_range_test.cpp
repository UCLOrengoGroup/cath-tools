/// \file
/// \brief The tuple_within_range test suite

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

#include "tuple_within_range.hpp"

#include <boost/test/unit_test.hpp>

using namespace ::cath::common;

using ::std::make_tuple;

BOOST_AUTO_TEST_SUITE(tuple_within_range_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert(   tuple_within_range( make_tuple(  0,  0,  0 ), make_tuple( 3, 7, 5 ) ) );
	static_assert(   tuple_within_range( make_tuple(  2,  6,  4 ), make_tuple( 3, 7, 5 ) ) );

	static_assert( ! tuple_within_range( make_tuple( -1,  0,  0 ), make_tuple( 3, 7, 5 ) ) );
	static_assert( ! tuple_within_range( make_tuple(  0, -1,  0 ), make_tuple( 3, 7, 5 ) ) );
	static_assert( ! tuple_within_range( make_tuple(  0,  0, -1 ), make_tuple( 3, 7, 5 ) ) );

	static_assert( ! tuple_within_range( make_tuple(  3,  6,  4 ), make_tuple( 3, 7, 5 ) ) );
	static_assert( ! tuple_within_range( make_tuple(  2,  7,  4 ), make_tuple( 3, 7, 5 ) ) );
	static_assert( ! tuple_within_range( make_tuple(  2,  6,  5 ), make_tuple( 3, 7, 5 ) ) );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
