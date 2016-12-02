/// \file
/// \brief The executable_options test suite

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
#include "executable_options.hpp"

#include <boost/test/unit_test.hpp>

using namespace cath::opts;

namespace cath {
	namespace test {

		struct executable_options_test_suite_fixture {
		protected:
			~executable_options_test_suite_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(executable_options_test_suite, cath::test::executable_options_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK( true );
}
BOOST_AUTO_TEST_SUITE_END()
