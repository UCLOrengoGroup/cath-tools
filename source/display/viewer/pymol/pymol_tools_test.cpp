/// \file
/// \brief The pymol_tools test suite

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

#include "common/size_t_literal.hpp"
#include "display/viewer/pymol/pymol_tools.hpp"
#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::common;

namespace cath {
	namespace test {

		/// \brief The pymol_tools_test_suite_fixture to assist in testing pymol_tools
		struct pymol_tools_test_suite_fixture : protected global_test_constants {
		protected:
			~pymol_tools_test_suite_fixture() noexcept = default;

			// Check that the asking for the y of either of the original x values gives back the original x value
			void check_pymol_size_ends(const size_t &arg_x_1,
			                           const double       &arg_y_1,
			                           const size_t &arg_x_2,
			                           const double       &arg_y_2
			                           ) {
				BOOST_CHECK_CLOSE( arg_y_1, pymol_tools::pymol_size(arg_x_1, arg_y_1, arg_x_2, arg_y_2, arg_x_1    ), ACCURACY_PERCENTAGE() );
				BOOST_CHECK_CLOSE( arg_y_2, pymol_tools::pymol_size(arg_x_1, arg_y_1, arg_x_2, arg_y_2, arg_x_2    ), ACCURACY_PERCENTAGE() );
				BOOST_CHECK_LT   ( 0.0,     pymol_tools::pymol_size(arg_x_1, arg_y_1, arg_x_2, arg_y_2, 10000000_z )                        );
			}
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(pymol_tools_test_suite, cath::test::pymol_tools_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(pymol_size) {
	check_pymol_size_ends( 1, 1.50, 100, 0.090 );
	check_pymol_size_ends( 1, 0.40, 100, 0.024 );
	check_pymol_size_ends( 1, 1.35, 100, 0.081 );
	check_pymol_size_ends( 1, 0.25, 100, 0.015 );
	check_pymol_size_ends( 1, 0.20, 100, 0.036 );
	check_pymol_size_ends( 1, 2.00, 100, 0.120 );
}

BOOST_AUTO_TEST_SUITE_END()

