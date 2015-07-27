/// \file
/// \brief The logger test suite

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

#include "logger.h"

#include <boost/test/unit_test.hpp>

namespace cath {
	namespace test {

		struct logger_fixture {
		protected:
			~logger_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(logger_test_suite, cath::test::logger_fixture)


/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL(1, 1);
}

BOOST_AUTO_TEST_SUITE_END()
