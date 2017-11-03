/// \file
/// \brief The coord_list class

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

#include <boost/test/auto_unit_test.hpp>

#include "structure/geometry/coord.hpp"
#include "structure/geometry/coord_list.hpp"

#include <vector>

using namespace cath;
using namespace cath::geom;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The coord_list_test_suite_fixture to assist in testing coord_list
		struct coord_list_test_suite_fixture {
		protected:
			~coord_list_test_suite_fixture() noexcept = default;

			coord_list coord_list_1 { coord_vec{
				{ 1.0, 2.0, 3.0 },
				{ 1.0, 2.0, 3.0 },
				{ 1.0, 2.0, 3.0 }
			} };
			coord_list coord_list_2 { coord_vec{
				{ 3.0, 2.0, 3.0 },
				{ 1.0, 4.0, 3.0 },
				{ 1.0, 2.0, 5.0 }
			} };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(coord_list_test_suite, cath::test::coord_list_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(rmsd) {
	BOOST_CHECK_EQUAL( 2.0, calc_rmsd( coord_list_1, coord_list_2 ) );
	BOOST_CHECK_EQUAL( 2.0, calc_rmsd( coord_list_2, coord_list_1 ) );
}

BOOST_AUTO_TEST_SUITE_END()
