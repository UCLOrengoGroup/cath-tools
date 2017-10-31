/// \file
/// \brief The residue test suite

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

#include "structure/protein/residue.hpp"

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "test/test_tools.hpp"

#include <ostream>

namespace cath { namespace common { class invalid_argument_exception; } }

using namespace cath;
using namespace cath::common;
using namespace cath::common::test;
using namespace cath::geom;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The residue_test_suite_fixture to assist in testing residue
		struct residue_test_suite_fixture {
		protected:
			~residue_test_suite_fixture() noexcept = default;

		public:
			residue get_residue_with_insert_without_ss() {
				return residue( make_residue_id( 'A', 43, 'A' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), residue::DEFAULT_PHI_PSI(), residue::DEFAULT_PHI_PSI(), 0 );
			}

			/// \brief A valine amino acid, used as the standard AA in the tests
			const amino_acid vaa { "VAL" };

			/// \brief A serine amino acid, used in the tests as an AA that's different from valine
			const amino_acid saa { "SER" };

			/// \brief An example residue for use in the tests
			///
			/// This isn't static because its construction depends on other statics,
			/// which might not be initialised yet
			const residue residue_with_insert_without_ss = { get_residue_with_insert_without_ss() };
		};

	}  // namespace test
}  // namespace cath

//
//const residue residue_test_suite_fixture::residue_with_insert_without_ss(
//	43, 'A', 'V', coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION, 0, 0, 0
//);

BOOST_FIXTURE_TEST_SUITE(residue_test_suite, cath::test::residue_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(invalid_insert_throws) {
	BOOST_CHECK_THROW(
		residue( make_residue_id( 'A', 43, '#' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), residue::DEFAULT_PHI_PSI(), residue::DEFAULT_PHI_PSI(), 0 ),
		invalid_argument_exception
	);
	BOOST_CHECK_THROW(
		residue( make_residue_id( 'A', 43, '!' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), residue::DEFAULT_PHI_PSI(), residue::DEFAULT_PHI_PSI(), 0 ),
		invalid_argument_exception
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(invalid_amino_acid_throws) {
	BOOST_CHECK_THROW(
		residue( make_residue_id( 'A', 43, 'A' ), amino_acid( '!' ), coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), residue::DEFAULT_PHI_PSI(), residue::DEFAULT_PHI_PSI(), 0 ),
		invalid_argument_exception
	);
	BOOST_CHECK_THROW(
		residue( make_residue_id( 'A', 43, 'A' ), amino_acid( 'a' ), coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), residue::DEFAULT_PHI_PSI(), residue::DEFAULT_PHI_PSI(), 0 ),
		invalid_argument_exception
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(output_empty_residue_left) {
	BOOST_CHECK_EQUAL("   0 0 0 0", ssap_legacy_alignment_left_side_gap_string());
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(output_empty_residue_right) {
	BOOST_CHECK_EQUAL("0 0 0    0", ssap_legacy_alignment_right_side_gap_string());
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(output_null_residue_left) {
	BOOST_CHECK_EQUAL("   0 0 0 X", ssap_legacy_alignment_left_side_string(residue::NULL_RESIDUE));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(output_null_residue_right) {
	BOOST_CHECK_EQUAL("X 0 0    0", ssap_legacy_alignment_right_side_string(residue::NULL_RESIDUE));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(output_residue_wiht_insert_without_ss_left) {
	BOOST_CHECK_EQUAL("  43 H A V", ssap_legacy_alignment_left_side_string(residue_with_insert_without_ss));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(output_residue_with_insert_without_ss_right) {
	BOOST_CHECK_EQUAL("V A H   43", ssap_legacy_alignment_right_side_string(residue_with_insert_without_ss));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(equality_works) {
	const auto     def_angle = residue::DEFAULT_PHI_PSI();
	const rotation other_rot = rotation::ROTATE_X_TO_Y_TO_Z_TO_X();
	const auto     svn_angle = make_angle_from_degrees<double>( 7.0 );
	const auto     egt_angle = make_angle_from_degrees<double>( 8.0 );
	const auto residues = {
		residue_with_insert_without_ss,
		residue(make_residue_id( 'A', 44, 'A' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43      ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), saa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), vaa, coord(2.2, 3.23, 4), coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), vaa, coord::ORIGIN_COORD, coord(2.2, 3.23, 4), 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 2, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::BETA_STRAND, rotation::IDENTITY_ROTATION(), def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, other_rot,                     def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), svn_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), def_angle, egt_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), vaa, coord::ORIGIN_COORD, coord::ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION(), def_angle, def_angle, 9)
	};
	check_equality_operators_on_diff_vals_range( residues );
}

/// \brief TODOCUMENT
/// rotates to:
///  * align the vector from N to C with the x-axis
///    (so that the N and C have equal y and z values and the C has greater x)
///  * align the plane between CA, N and C with the x-z plane
///    (so that the CA's y value equals those of the N and C and its z value is greater the N and C have equal y and z values and the C has greater x)
BOOST_AUTO_TEST_CASE(construct_residue_frame_works) {
	BOOST_CHECK_EQUAL(
		rotation::IDENTITY_ROTATION(),
		construct_residue_frame(
			coord(0, 0, 0),
			coord(1, 0, 1),
			coord(2, 0, 0)
		)
	);
	BOOST_CHECK_EQUAL(
		rotation::ROTATE_X_TO_Y_TO_Z_TO_X(),
		construct_residue_frame(
			coord(0, 0, 0),
			coord(0, 1, 1),
			coord(0, 0, 2)
		)
	);
	BOOST_CHECK_EQUAL(
		rotation::ROTATE_X_TO_Z_TO_Y_TO_X(),
		construct_residue_frame(
			coord(0, 0, 0),
			coord(1, 1, 0),
			coord(0, 2, 0)
		)
	);
	BOOST_CHECK_EQUAL(
		rotation::ROTATE_X_TO_Z_TO_Y_TO_X(),
		construct_residue_frame(
			coord(5, 5, 5),
			coord(6, 6, 5),
			coord(5, 7, 5)
		)
	);
}

BOOST_AUTO_TEST_SUITE_END()
