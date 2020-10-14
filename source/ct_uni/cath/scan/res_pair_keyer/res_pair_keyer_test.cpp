/// \file
/// \brief The res_pair_keyer test suite

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
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

#include <boost/test/unit_test.hpp>

#include "cath/scan/res_pair_keyer/res_pair_keyer.hpp"

// #include "cath/test/global_test_constants.hpp"

using namespace ::cath::scan::detail;
//using namespace ::std;

namespace cath {
	namespace test {

		/// \brief The res_pair_keyer_test_suite_fixture to assist in testing res_pair_keyer
		struct res_pair_keyer_test_suite_fixture {
		protected:
			~res_pair_keyer_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(res_pair_keyer_test_suite, cath::test::res_pair_keyer_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
