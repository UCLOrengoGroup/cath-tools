/// \file
/// \brief The tuple_lattice_index test suite

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

#include "tuple_lattice_index.hpp"

#include <boost/test/unit_test.hpp>

// #include <array>

using namespace cath::common;

// using std::array;
using std::make_tuple;
// using std::tuple;

BOOST_AUTO_TEST_SUITE(tuple_lattice_index_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert( tuple_lattice_index( make_tuple(          2 ), make_tuple(          3 ) ) ==  2, "tuple_lattice_index() should give this result" );
	static_assert( tuple_lattice_index( make_tuple(          5 ), make_tuple(          7 ) ) ==  5, "tuple_lattice_index() should give this result" );
	static_assert( tuple_lattice_index( make_tuple(         11 ), make_tuple(         13 ) ) == 11, "tuple_lattice_index() should give this result" );

	static_assert( tuple_lattice_index( make_tuple(      2,  5 ), make_tuple(      3,  7 ) ) == 19, "tuple_lattice_index() should give this result" );
	static_assert( tuple_lattice_index( make_tuple(      2, 11 ), make_tuple(      3, 13 ) ) == 37, "tuple_lattice_index() should give this result" );
	static_assert( tuple_lattice_index( make_tuple(      5, 11 ), make_tuple(      7, 13 ) ) == 76, "tuple_lattice_index() should give this result" );
	static_assert( tuple_lattice_index( make_tuple(      5,  2 ), make_tuple(      7,  3 ) ) == 17, "tuple_lattice_index() should give this result" );
	static_assert( tuple_lattice_index( make_tuple(     11,  2 ), make_tuple(     13,  3 ) ) == 35, "tuple_lattice_index() should give this result" );
	static_assert( tuple_lattice_index( make_tuple(     11,  5 ), make_tuple(     13,  7 ) ) == 82, "tuple_lattice_index() should give this result" );

	static_assert( tuple_lattice_index( make_tuple(  2,  5, 11 ), make_tuple(  3,  7, 13 ) ) == 258, "tuple_lattice_index() should give this result" );
	static_assert( tuple_lattice_index( make_tuple( 11,  5,  2 ), make_tuple( 13,  7,  3 ) ) == 248, "tuple_lattice_index() should give this result" );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
