/// \file
/// \brief The superposition test suite

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

#include <boost/test/output_test_stream.hpp>

#include "common/test_tools.h"
#include "test/global_test_constants.h"
#include "structure/geometry/coord.h"
#include "structure/geometry/coord_list.h"
#include "structure/structure_type_aliases.h"
#include "superposition/superposition.h"
#include "superposition/io/superposition_io.h"

//#include <iostream> // *** TEMPORARY ***
#include <vector>

using namespace boost::test_tools;
using namespace cath;
using namespace cath::common::test;
using namespace cath::geom;
using namespace cath::sup;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The superposition_test_suite_fixture to assist in testing superposition
		struct superposition_test_suite_fixture : public global_test_constants {
		protected:
			~superposition_test_suite_fixture() noexcept = default;

		public:
			const coord_list coord_list1 { coord_vec{
				{ 10.306, 134.301, 48.038 },
				{ 13.342, 136.603, 47.735 },
				{ 16.498, 135.180, 46.225 },
				{ 19.855, 136.695, 45.235 },
				{ 21.726, 135.373, 42.223 }
			} };

			const coord_list coord_list2 { coord_vec{
				{ 21.222, 18.395, 14.116 },
				{ 24.152, 16.183, 15.172 },
				{ 23.817, 12.426, 14.864 },
				{ 25.966,  9.358, 15.175 },
				{ 25.362,  6.469, 12.836 }
			} };

			const coord_list coord_list2_variation { coord_vec{
				{ 21.222, 18.395, 14.116 },
				{ 24.152, 16.183, 15.172 },
				{ 23.817, 12.426, 14.864 },
				{ 25.966,  9.358, 15.175 },
				{ 15.000, 15.000, 15.000 }
			} };

			const double rmsd_between_1_and_2 = { 0.12978869963736103 };
		};
	}
}

BOOST_FIXTURE_TEST_SUITE(superposition_test_suite, cath::test::superposition_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(equality) {
	const auto superpositions = {
		create_pairwise_superposition( coord_list1, coord_list2           ),
		create_pairwise_superposition( coord_list1, coord_list2_variation )
	};
	check_equality_operators_on_diff_vals_range( superpositions );
}

/// \brief Check that things work as expected if specifying:
///         - a different entry to be used as the base
///         - a translation to be applied to the base structure
///         - a rotation    to be applied to the base structure
BOOST_AUTO_TEST_CASE(vary_base_rot_and_trans) {
	const superposition test_sup_first         = create_pairwise_superposition( coord_list1, coord_list2, true  );
	const coord         sup_first_translation  = test_sup_first.get_translation_of_index(  superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION );
	const rotation      sup_first_rotation     = test_sup_first.get_rotation_of_index(     superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION );

	const superposition test_sup_second        = create_pairwise_superposition( coord_list1, coord_list2, false );
	const coord         sup_second_translation = test_sup_second.get_translation_of_index( superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION  );
	const rotation      sup_second_rotation    = test_sup_second.get_rotation_of_index(    superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION  );

	const superposition test_sup_first_rev     = create_pairwise_superposition( coord_list1, coord_list2, true,  sup_second_translation, sup_second_rotation );
	const superposition test_sup_second_rev    = create_pairwise_superposition( coord_list1, coord_list2, false, sup_first_translation,  sup_first_rotation  );

	BOOST_CHECK(  are_close( test_sup_first,  test_sup_second_rev ) );
	BOOST_CHECK(  are_close( test_sup_second, test_sup_first_rev  ) );
	BOOST_CHECK( !are_close( test_sup_first,  test_sup_second     ) );
}

/// \brief Check that superposition can handle one coord_list being a perfect 180 degree rotation of the other
///
/// \todo Have this testcase generate the flipped versions rather than hard-coding them in here
BOOST_AUTO_TEST_CASE(axis_pair_negation_y_and_z) {
	const coord_list coord_list1_flipped_around_x{ {
		{ 10.306, -134.301, -48.038 },
		{ 13.342, -136.603, -47.735 },
		{ 16.498, -135.180, -46.225 },
		{ 19.855, -136.695, -45.235 },
		{ 21.726, -135.373, -42.223 }
	} };
	const coord_list coord_list1_flipped_around_y{ {
		{ -10.306,  134.301, -48.038 },
		{ -13.342,  136.603, -47.735 },
		{ -16.498,  135.180, -46.225 },
		{ -19.855,  136.695, -45.235 },
		{ -21.726,  135.373, -42.223 }
	} };
	const coord_list coord_list1_flipped_around_z{ {
		{ -10.306, -134.301,  48.038 },
		{ -13.342, -136.603,  47.735 },
		{ -16.498, -135.180,  46.225 },
		{ -19.855, -136.695,  45.235 },
		{ -21.726, -135.373,  42.223 }
	} };
	const double RMSD_ACCURACY(0.00000001);
	BOOST_CHECK_SMALL(calc_pairwise_superposition_rmsd(coord_list1, coord_list1_flipped_around_x), RMSD_ACCURACY );
	BOOST_CHECK_SMALL(calc_pairwise_superposition_rmsd(coord_list1, coord_list1_flipped_around_y), RMSD_ACCURACY );
	BOOST_CHECK_SMALL(calc_pairwise_superposition_rmsd(coord_list1, coord_list1_flipped_around_z), RMSD_ACCURACY );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(rmsd) {
	BOOST_CHECK_CLOSE(
		rmsd_between_1_and_2,
		calc_pairwise_superposition_rmsd(coord_list1, coord_list2),
		ACCURACY_PERCENTAGE()
	);
}

BOOST_AUTO_TEST_CASE(to_json_string_works_for_example_sup) {
	BOOST_CHECK_EQUAL(
		to_json_string( create_pairwise_superposition( coord_list1, coord_list2 ), false ),
		R"({"transformations":[{"translation":)"
		R"({"x":"0","y":"0","z":"0"},)"
		R"("rotation":)"
		R"([["1","0","0"],)"
		R"(["0","1","0"],)"
		R"(["0","0","1"])"
		R"(]},{"translation":)"
		R"({"x":"104.47726177086807","y":"20.66883749985081","z":"41.523834652502032"},)"
		R"("rotation":)"
		R"([["0.20769400083047135","-0.9111116114406772","0.35600397963647046"],)"
		R"(["0.96878664188870112","0.24194236224656662","0.0540031096194028"],)"
		R"(["-0.13533530403056787","0.33367577803689003","0.93292273561878114"])"
		R"(]}]})"
		"\n"
	);
}


BOOST_AUTO_TEST_SUITE_END()

