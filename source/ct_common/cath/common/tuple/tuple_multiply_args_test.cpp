/// \file
/// \brief The tuple_multiply_args test suite

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

#include "tuple_multiply_args.hpp"

#include <boost/test/unit_test.hpp>

using namespace ::cath::common;

using ::std::make_tuple;

BOOST_AUTO_TEST_SUITE(tuple_multiply_args_test_suite)

BOOST_AUTO_TEST_CASE(basic) {

	static_assert( tuple_multiply_args( make_tuple( 3,   7,   2   ) ) == 42 );
	static_assert( tuple_multiply_args( make_tuple( 3.0, 7.0, 2.0 ) ) == 42.0 );

	// static_assert( tuple_multiply_args( ) == 0, "" ); // This should fail to compile with something like: `no matching function for call to object of type 'tuple_multiply_args_fn'`

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
