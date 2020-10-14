/// \file
/// \brief The make_tuple_with_skips test suite

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

#include "make_tuple_with_skips.hpp"

#include <boost/test/unit_test.hpp>

using namespace ::cath::common;

using ::std::make_tuple;

BOOST_AUTO_TEST_SUITE(make_tuple_with_skips_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert(
		make_tuple_with_skips( 5, tpl_elmnt_skip_t{}, 4u )
		==
		make_tuple           ( 5,                     4u ),
		""
	);

	constexpr tpl_elmnt_skip_t my_skipper{};
	static_assert(
		make_tuple_with_skips( 5, my_skipper,         4u )
		==
		make_tuple           ( 5,                     4u ),
		""
	);

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
