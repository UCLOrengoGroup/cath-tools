/// \file
/// \brief The temp_check_offset_1 test suite

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

#include "temp_check_offset_1.hpp"

#include <boost/test/unit_test.hpp>

#include "cath/common/config.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"

using namespace ::cath;
using namespace ::cath::common;

BOOST_AUTO_TEST_SUITE(temp_check_offset_1_test_suite)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE( zero_throws ) {
	if constexpr ( IS_IN_DEBUG_MODE ) {
		BOOST_CHECK_THROW( check_offset_1( 0 ), invalid_argument_exception );
	} else {
		BOOST_CHECK_NO_THROW_DIAG( check_offset_1( 0 ) );
	}
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(one_does_not_throw) {
	BOOST_CHECK_NO_THROW_DIAG( check_offset_1(1) );
}

BOOST_AUTO_TEST_SUITE_END()
