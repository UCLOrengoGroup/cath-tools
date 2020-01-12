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

#include "structure/geometry/coord.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/structure_type_aliases.hpp"
#include "superposition/io/superposition_io.hpp"
#include "superposition/superposition.hpp"
#include "test/global_test_constants.hpp"
#include "test/test_tools.hpp"

//#include <iostream> // *** TEMPORARY ***
#include <vector>

using namespace cath;
using namespace cath::common::test;
using namespace cath::geom;
using namespace cath::sup;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The superposition_test_suite_fixture to assist in testing superposition
		struct superposition_test_suite_fixture : protected global_test_constants {
		protected:
			~superposition_test_suite_fixture() noexcept = default;

		public:
			const coord_list coord_list_1 { coord_vec{
				{ 10.306, 134.301, 48.038 },
				{ 13.342, 136.603, 47.735 },
				{ 16.498, 135.180, 46.225 },
				{ 19.855, 136.695, 45.235 },
				{ 21.726, 135.373, 42.223 }
			} };

			const coord_list coord_list_2 { coord_vec{
				{ 21.222, 18.395, 14.116 },
				{ 24.152, 16.183, 15.172 },
				{ 23.817, 12.426, 14.864 },
				{ 25.966,  9.358, 15.175 },
				{ 25.362,  6.469, 12.836 }
			} };

			const coord_list coord_list_2_variation { coord_vec{
				{ 21.222, 18.395, 14.116 },
				{ 24.152, 16.183, 15.172 },
				{ 23.817, 12.426, 14.864 },
				{ 25.966,  9.358, 15.175 },
				{ 15.000, 15.000, 15.000 }
			} };

			const double rmsd_between_1_and_2 = { 0.12978869963736103 };

			/// \brief Check that superposition returned by post_translate_and_rotate() has the
			///        same effect on 7 key points as applying the parts separately
			void check_post_translate_and_rotate(const coord    &prm_orig_supn_transltn, ///< The translation part of the original superposition to apply first
			                                     const rotation &prm_orig_supn_rottn,    ///< The translation part of the original superposition to apply first
			                                     const coord    &prm_translation,        ///< The translation to apply after the superposition
			                                     const rotation &prm_rotation            ///< The rotation to apply last, after the translation
			                                     ) {
				constexpr size_t IDX = 0;
				const superposition orig_supn{ { prm_orig_supn_transltn }, { prm_orig_supn_rottn } };
				const superposition ptar_supn = post_translate_and_rotate_copy( orig_supn, prm_translation, prm_rotation );
				for (const coord &x : {  ORIGIN_COORD,
				                         UNIT_X_COORD,  UNIT_Y_COORD,  UNIT_Z_COORD,
				                        -UNIT_X_COORD, -UNIT_Y_COORD, -UNIT_Z_COORD } ) {
					BOOST_CHECK_EQUAL(
						transform_copy( ptar_supn, IDX, x ),
						rotate_copy(
							prm_rotation,
							transform_copy( orig_supn, IDX, x ) + prm_translation
						)
					);
				}
			}
		};
	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(superposition_test_suite, cath::test::superposition_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(equality) {
	const auto superpositions = {
		create_pairwise_superposition( coord_list_1, coord_list_2           ),
		create_pairwise_superposition( coord_list_1, coord_list_2_variation )
	};
	check_equality_operators_on_diff_vals_range( superpositions );
}

/// \brief Check that things work as expected if specifying:
///         - a different entry to be used as the base
///         - a translation to be applied to the base structure
///         - a rotation    to be applied to the base structure
BOOST_AUTO_TEST_CASE(vary_base_rot_and_trans) {
	const superposition test_sup_first         = create_pairwise_superposition( coord_list_1, coord_list_2, true  );
	const coord         sup_first_translation  = test_sup_first.get_translation_of_index(  superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION );
	const rotation      sup_first_rotation     = test_sup_first.get_rotation_of_index(     superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION );

	const superposition test_sup_second        = create_pairwise_superposition( coord_list_1, coord_list_2, false );
	const coord         sup_second_translation = test_sup_second.get_translation_of_index( superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION  );
	const rotation      sup_second_rotation    = test_sup_second.get_rotation_of_index(    superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION  );

	const superposition test_sup_first_rev     = create_pairwise_superposition( coord_list_1, coord_list_2, true,  sup_second_translation, sup_second_rotation );
	const superposition test_sup_second_rev    = create_pairwise_superposition( coord_list_1, coord_list_2, false, sup_first_translation,  sup_first_rotation  );

	BOOST_CHECK(  are_close( test_sup_first,  test_sup_second_rev ) );
	BOOST_CHECK(  are_close( test_sup_second, test_sup_first_rev  ) );
	BOOST_CHECK( !are_close( test_sup_first,  test_sup_second     ) );
}

/// \brief Check that superposition can handle one coord_list being a perfect 180 degree rotation of the other
///
/// \todo Have this testcase generate the flipped versions rather than hard-coding them in here
BOOST_AUTO_TEST_CASE(axis_pair_negation_y_and_z) {
	const coord_list coord_list_1_flipped_around_x{ {
		{ 10.306, -134.301, -48.038 },
		{ 13.342, -136.603, -47.735 },
		{ 16.498, -135.180, -46.225 },
		{ 19.855, -136.695, -45.235 },
		{ 21.726, -135.373, -42.223 }
	} };
	const coord_list coord_list_1_flipped_around_y{ {
		{ -10.306,  134.301, -48.038 },
		{ -13.342,  136.603, -47.735 },
		{ -16.498,  135.180, -46.225 },
		{ -19.855,  136.695, -45.235 },
		{ -21.726,  135.373, -42.223 }
	} };
	const coord_list coord_list_1_flipped_around_z{ {
		{ -10.306, -134.301,  48.038 },
		{ -13.342, -136.603,  47.735 },
		{ -16.498, -135.180,  46.225 },
		{ -19.855, -136.695,  45.235 },
		{ -21.726, -135.373,  42.223 }
	} };
	const double RMSD_ACCURACY(0.00000001);
	BOOST_CHECK_SMALL(calc_pairwise_superposition_rmsd(coord_list_1, coord_list_1_flipped_around_x), RMSD_ACCURACY );
	BOOST_CHECK_SMALL(calc_pairwise_superposition_rmsd(coord_list_1, coord_list_1_flipped_around_y), RMSD_ACCURACY );
	BOOST_CHECK_SMALL(calc_pairwise_superposition_rmsd(coord_list_1, coord_list_1_flipped_around_z), RMSD_ACCURACY );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(rmsd) {
	BOOST_CHECK_CLOSE(
		rmsd_between_1_and_2,
		calc_pairwise_superposition_rmsd(coord_list_1, coord_list_2),
		ACCURACY_PERCENTAGE()
	);
}

BOOST_AUTO_TEST_CASE(post_translate_and_rotate_work) {
	const auto rotations = {
		rotation::ROTATE_X_TO_Y_TO_Z_TO_X(),
		rotation::ROTATE_X_TO_Z_TO_Y_TO_X(),
	};
	const auto translations = {
		coord{ 1.0, 2.0, 3.0 },
		coord{ 3.0, 1.0, 2.0 },
		coord{ 2.0, 3.0, 1.0 },
	};
	for (const coord &orig_trans : translations) {
		for (const rotation &orig_rotn : rotations) {
			for (const coord &post_trans : translations) {
				for (const rotation &post_rotn : rotations) {
					check_post_translate_and_rotate( orig_trans, orig_rotn, post_trans, post_rotn );
				}
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()

