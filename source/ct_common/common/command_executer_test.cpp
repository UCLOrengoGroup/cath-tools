/// \file
/// \brief The command_executer test suite

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

namespace cath {
	namespace test {

		/// \brief The command_executer_test_suite_fixture to assist in testing command_executer
		struct command_executer_test_suite_fixture {
		protected:
			~command_executer_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(command_executer_test_suite, cath::test::command_executer_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK( true );
}
BOOST_AUTO_TEST_SUITE_END()

