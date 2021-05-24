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

#include "cath/structure/protein/residue.hpp"

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "cath/test/test_tools.hpp"

#include <ostream>

// clang-format off
namespace cath::common { class invalid_argument_exception; }
// clang-format on

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::common::test;
using namespace ::cath::geom;
using namespace ::std;

/// \brief A valine amino acid, used as the standard AA in the tests
constexpr amino_acid VAA{ "VAL" };

/// \brief A serine amino acid, used in the tests as an AA that's different from valine
constexpr amino_acid SAA{ "SER" };

/// \brief An example residue for use in the tests
constexpr residue RESIDUE_WITH_INSERT_WITHOUT_SS{ make_residue_id( 'A', 43, 'A' ),
                                                  VAA,
                                                  ORIGIN_COORD,
                                                  ORIGIN_COORD,
                                                  3,
                                                  sec_struc_type::ALPHA_HELIX,
                                                  IDENTITY_ROTATION,
                                                  residue::DEFAULT_PHI_PSI,
                                                  residue::DEFAULT_PHI_PSI,
                                                  0 };

//
//const residue residue_test_suite_fixture::RESIDUE_WITH_INSERT_WITHOUT_SS(
//	43, 'A', 'V', ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, rotation::IDENTITY_ROTATION, 0, 0, 0
//);

BOOST_AUTO_TEST_SUITE(residue_test_suite)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(invalid_insert_throws) {
	BOOST_CHECK_THROW(
		residue( make_residue_id( 'A', 43, '#' ), VAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, residue::DEFAULT_PHI_PSI, residue::DEFAULT_PHI_PSI, 0 ),
		invalid_argument_exception
	);
	BOOST_CHECK_THROW(
		residue( make_residue_id( 'A', 43, '!' ), VAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, residue::DEFAULT_PHI_PSI, residue::DEFAULT_PHI_PSI, 0 ),
		invalid_argument_exception
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(invalid_amino_acid_throws) {
	BOOST_CHECK_THROW(
		residue( make_residue_id( 'A', 43, 'A' ), amino_acid( '!' ), ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, residue::DEFAULT_PHI_PSI, residue::DEFAULT_PHI_PSI, 0 ),
		invalid_argument_exception
	);
	BOOST_CHECK_THROW(
		residue( make_residue_id( 'A', 43, 'A' ), amino_acid( 'a' ), ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, residue::DEFAULT_PHI_PSI, residue::DEFAULT_PHI_PSI, 0 ),
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
	BOOST_CHECK_EQUAL("   0 0 0 X", ssap_legacy_alignment_left_side_string(NULL_RESIDUE));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(output_null_residue_right) {
	BOOST_CHECK_EQUAL("X 0 0    0", ssap_legacy_alignment_right_side_string(NULL_RESIDUE));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(output_residue_wiht_insert_without_ss_left) {
	BOOST_CHECK_EQUAL("  43 H A V", ssap_legacy_alignment_left_side_string(RESIDUE_WITH_INSERT_WITHOUT_SS));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(output_RESIDUE_WITH_INSERT_WITHOUT_SS_right) {
	BOOST_CHECK_EQUAL("V A H   43", ssap_legacy_alignment_right_side_string(RESIDUE_WITH_INSERT_WITHOUT_SS));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(equality_works) {
	const auto     def_angle = residue::DEFAULT_PHI_PSI;
	const rotation other_rot = ROTATE_X_TO_Y_TO_Z_TO_X;
	const auto     svn_angle = make_angle_from_degrees<double>( 7.0 );
	const auto     egt_angle = make_angle_from_degrees<double>( 8.0 );
	const auto residues = {
		RESIDUE_WITH_INSERT_WITHOUT_SS,
		residue(make_residue_id( 'A', 44, 'A' ), VAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43      ), VAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), SAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), VAA, coord(2.2, 3.23, 4), ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), VAA, ORIGIN_COORD, coord(2.2, 3.23, 4), 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), VAA, ORIGIN_COORD, ORIGIN_COORD, 2, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), VAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::BETA_STRAND, IDENTITY_ROTATION, def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), VAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, other_rot,                     def_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), VAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, svn_angle, def_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), VAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, def_angle, egt_angle, 0),
		residue(make_residue_id( 'A', 43, 'A' ), VAA, ORIGIN_COORD, ORIGIN_COORD, 3, sec_struc_type::ALPHA_HELIX, IDENTITY_ROTATION, def_angle, def_angle, 9)
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
		IDENTITY_ROTATION,
		construct_residue_frame(
			coord(0, 0, 0),
			coord(1, 0, 1),
			coord(2, 0, 0)
		)
	);
	BOOST_CHECK_EQUAL(
		ROTATE_X_TO_Y_TO_Z_TO_X,
		construct_residue_frame(
			coord(0, 0, 0),
			coord(0, 1, 1),
			coord(0, 0, 2)
		)
	);
	BOOST_CHECK_EQUAL(
		ROTATE_X_TO_Z_TO_Y_TO_X,
		construct_residue_frame(
			coord(0, 0, 0),
			coord(1, 1, 0),
			coord(0, 2, 0)
		)
	);
	BOOST_CHECK_EQUAL(
		ROTATE_X_TO_Z_TO_Y_TO_X,
		construct_residue_frame(
			coord(5, 5, 5),
			coord(6, 6, 5),
			coord(5, 7, 5)
		)
	);
}

BOOST_AUTO_TEST_SUITE_END()
