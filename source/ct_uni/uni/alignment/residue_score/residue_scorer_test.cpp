/// \file
/// \brief The residue_scorer test suite

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

#include "alignment/residue_score/residue_scorer.hpp"
#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::align;

namespace cath {
	namespace test {

		/// \brief The residue_scorer_test_suite_fixture to assist in testing residue_scorer
		struct residue_scorer_test_suite_fixture: protected global_test_constants {
		protected:
			~residue_scorer_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(residue_scorer_test_suite, cath::test::residue_scorer_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
