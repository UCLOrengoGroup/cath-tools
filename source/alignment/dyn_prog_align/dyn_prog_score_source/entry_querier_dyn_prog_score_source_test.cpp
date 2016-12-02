/// \file
/// \brief The entry_querier_dyn_prog_score_source test suite

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

#include "entry_querier_dyn_prog_score_source.hpp"

#include <boost/test/unit_test.hpp>

namespace cath {
	namespace test {

		struct entry_querier_dyn_prog_score_source_fixture {
		protected:
			~entry_querier_dyn_prog_score_source_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(entry_querier_dyn_prog_score_source_test_suite, cath::test::entry_querier_dyn_prog_score_source_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
