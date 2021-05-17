/// \file
/// \brief The make_optional_if test suite

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

#include "make_optional_if.hpp"

#include <boost/test/unit_test.hpp>

#include "cath/test/boost_test_print_type.hpp"

using namespace ::cath::common;

using ::std::make_optional;
using ::std::nullopt;
using ::std::optional;

BOOST_AUTO_TEST_SUITE(make_optional_if_test_suite)

BOOST_AUTO_TEST_CASE(make_optional_if_works) {
	BOOST_TEST( make_optional_if( true,  1 ) == 1       );
	BOOST_TEST( make_optional_if( false, 1 ) == nullopt );
}

BOOST_AUTO_TEST_CASE(make_optional_if_fn_works) {
	static_assert( make_optional_if_fn( true,  [] {                      return 1; } ) == 1       );
	static_assert( make_optional_if_fn( false, [] { BOOST_CHECK( false); return 1; } ) == nullopt );
}

BOOST_AUTO_TEST_CASE( if_then_optional_works ) {
	int a = 0;
	BOOST_TEST( if_then_optional( true, ++a ) == 1 );
	BOOST_TEST( a == 1 );
	BOOST_TEST( if_then_optional( false, ++a ) == nullopt );
	BOOST_TEST( a == 1 );
}

BOOST_AUTO_TEST_SUITE_END()
