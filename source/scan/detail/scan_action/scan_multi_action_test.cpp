/// \file
/// \brief The scan_multi_action test suite

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
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

#include <boost/test/auto_unit_test.hpp>

#include "scan/detail/scan_action/scan_multi_action.h"

// #include "test/global_test_constants.h"

using namespace cath::scan::detail;
//using namespace std;

namespace cath {
	namespace test {

		/// \brief The scan_multi_action_test_suite_fixture to assist in testing scan_multi_action
		struct scan_multi_action_test_suite_fixture {
		protected:
			~scan_multi_action_test_suite_fixture() noexcept = default;
		};

	}
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(scan_multi_action_test_suite, cath::test::scan_multi_action_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL(0, 0);
}

BOOST_AUTO_TEST_SUITE_END()
