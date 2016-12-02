/// \file
/// \brief The variadic_and test suite

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

#include "variadic_and.h"

#include <boost/test/auto_unit_test.hpp>

using namespace cath::common;

BOOST_AUTO_TEST_SUITE(variadic_and_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert( variadic_and( true  ) == true,  "" );
	static_assert( variadic_and( false ) == false, "" );

	static_assert( variadic_and( true,  true  ) == true,  "" );
	static_assert( variadic_and( true,  false ) == false, "" );
	static_assert( variadic_and( false, true  ) == false, "" );
	static_assert( variadic_and( false, false ) == false, "" );

	static_assert( variadic_and( true,  true,  true  ) == true,  "" );
	static_assert( variadic_and( true,  true,  false ) == false, "" );
	static_assert( variadic_and( true,  false, true  ) == false, "" );
	static_assert( variadic_and( true,  false, false ) == false, "" );
	static_assert( variadic_and( false, true,  true  ) == false, "" );
	static_assert( variadic_and( false, true,  false ) == false, "" );
	static_assert( variadic_and( false, false, true  ) == false, "" );
	static_assert( variadic_and( false, false, false ) == false, "" );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
