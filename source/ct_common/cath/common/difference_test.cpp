/// \file
/// \brief The difference test suite

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

#include <boost/test/unit_test.hpp>

#include "difference.hpp"

#include "cath/common/size_t_literal.hpp"

#include <cmath>

using namespace ::cath::common;
using namespace ::std;

/// \brief TODOCUMENT
BOOST_AUTO_TEST_SUITE(difference_test_suite)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	const size_t five( 5_z );
	const size_t nine( 9_z );
	BOOST_CHECK_EQUAL( difference( five,  nine ), 4_z );
	BOOST_CHECK_EQUAL( difference( nine,  five ), 4_z );
}

BOOST_AUTO_TEST_SUITE_END()
