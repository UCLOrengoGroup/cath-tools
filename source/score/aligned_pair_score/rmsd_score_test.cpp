/// \file
/// \brief The rmsd_score test suite

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

#include "rmsd_score.h"

#include <boost/test/unit_test.hpp>

#include "alignment/common_atom_selection_policy/common_atom_selection_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.h"

namespace cath {
	namespace test {

		struct rmsd_score_fixture {
		protected:
			~rmsd_score_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(rmsd_score_test_suite, cath::test::rmsd_score_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( 0, 0 );
}

BOOST_AUTO_TEST_SUITE_END()