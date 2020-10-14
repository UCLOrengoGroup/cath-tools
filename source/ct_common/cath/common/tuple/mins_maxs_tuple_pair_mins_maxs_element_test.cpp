/// \file
/// \brief The mins_maxs_tuple_pair_mins_maxs_element test suite

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

#include "mins_maxs_tuple_pair_mins_maxs_element.hpp"

#include <boost/test/unit_test.hpp>

#include <array>

using namespace cath::common;

using std::array;
using std::make_tuple;
using std::pair;
using std::tuple;

BOOST_AUTO_TEST_SUITE(mins_maxs_tuple_pair_mins_maxs_element_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	
	constexpr array< pair< tuple<int, int, int>, tuple<int, int, int> >, 5> bob{ {
		make_pair( make_tuple( -29, -13,  72 ), make_tuple( -29, -13,  72 ) ),
		make_pair( make_tuple(  36, -78, -18 ), make_tuple(  36, -78, -18 ) ),
		make_pair( make_tuple( -45,  75, -60 ), make_tuple( -45,  75, -60 ) ),
		make_pair( make_tuple(  88, -47, -55 ), make_tuple(  88, -47, -55 ) ),
		make_pair( make_tuple(  68, -38, -58 ), make_tuple(  68, -38, -58 ) ),
	} };
	BOOST_CHECK_EQUAL( std::get<0>( mins_maxs_tuple_pair_mins_maxs_element( bob ).first  ), -45 );
	BOOST_CHECK_EQUAL( std::get<1>( mins_maxs_tuple_pair_mins_maxs_element( bob ).first  ), -78 );
	BOOST_CHECK_EQUAL( std::get<2>( mins_maxs_tuple_pair_mins_maxs_element( bob ).first  ), -60 );
	BOOST_CHECK_EQUAL( std::get<0>( mins_maxs_tuple_pair_mins_maxs_element( bob ).second ),  88 );
	BOOST_CHECK_EQUAL( std::get<1>( mins_maxs_tuple_pair_mins_maxs_element( bob ).second ),  75 );
	BOOST_CHECK_EQUAL( std::get<2>( mins_maxs_tuple_pair_mins_maxs_element( bob ).second ),  72 );
}

BOOST_AUTO_TEST_SUITE_END()
