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

#include <boost/optional/optional_io.hpp>
#include <boost/test/unit_test.hpp>

using namespace cath::common;

using boost::none;
using boost::optional;

BOOST_AUTO_TEST_SUITE(make_optional_if_test_suite)

BOOST_AUTO_TEST_CASE(make_optional_if_works) {
	BOOST_CHECK_EQUAL( make_optional_if( true,  1 ), optional<int>{ 1    } );
	BOOST_CHECK_EQUAL( make_optional_if( false, 1 ), optional<int>{ none } );
}

BOOST_AUTO_TEST_CASE(make_optional_if_fn_works) {
	BOOST_CHECK_EQUAL( make_optional_if_fn( true,  [] {                      return 1; } ), optional<int>{ 1    } );
	BOOST_CHECK_EQUAL( make_optional_if_fn( false, [] { BOOST_CHECK( false); return 1; } ), optional<int>{ none } );
}

BOOST_AUTO_TEST_SUITE_END()
