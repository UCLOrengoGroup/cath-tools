/// \file
/// \brief The sec_play test suite

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

#include "sec_play.hpp"

#include "common/boost_addenda/range/front.hpp"
#include "file/dssp_wolf/dssp_file.hpp"
#include "file/dssp_wolf/dssp_file_io.hpp"
#include "file/pdb/pdb.hpp"
#include "file/sec/sec_file.hpp"
#include "file/sec/sec_file_io.hpp"
#include "file/sec/sec_file_record.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sec;

BOOST_FIXTURE_TEST_SUITE(sec_play_test_suite, global_test_constants)

// Starts and stops on 1hdoA
// {   5,   9 },
// {  14,  25 },
// {  29,  34 },
// {  49,  52 },
// {  58,  65 },
// {  70,  73 },
// {  86, 101 },
// { 105, 109 },
// { 126, 141 },
// { 145, 149 },
// { 152, 156 },
// { 163, 167 },
// { 174, 178 },
// { 178, 187 },
// { 198, 202 },

BOOST_AUTO_TEST_CASE(get_correct_for_example_b_1hdoA00_first_strand) {
	const protein the_protein = protein_from_dssp_and_pdb(
		read_dssp_file( EXAMPLE_B_DSSP_FILENAME() ),
		read_pdb_file ( EXAMPLE_B_PDB_FILENAME()  )
	);
	const sec_file_record got_first_sec      = calculate_sec_record( the_protein, 4, 8        );
	const sec_file_record expected_first_sec = front( read_sec     ( EXAMPLE_B_SEC_FILENAME() ) );
	BOOST_CHECK_EQUAL( round_like_sec_file_copy( got_first_sec ), expected_first_sec );
}

BOOST_AUTO_TEST_CASE(get_correct_for_example_b_1hdoA00_later_helix) {
	const protein the_protein = protein_from_dssp_and_pdb(
		read_dssp_file( EXAMPLE_B_DSSP_FILENAME() ),
		read_pdb_file ( EXAMPLE_B_PDB_FILENAME()  )
	);
	const sec_file_record got_sec      = calculate_sec_record( the_protein, 57, 64 );
	const sec_file_record expected_sec = *next( common::cbegin( read_sec( EXAMPLE_B_SEC_FILENAME() ) ), 4 );
	BOOST_CHECK_EQUAL( round_like_sec_file_copy( got_sec ), expected_sec );
}

BOOST_AUTO_TEST_CASE(prosec_axis_point_gets_beta_strand_example_correct) {
	static_assert(
		prosec_axis_point_of_residue_triple(
			coord{ -12.743,  -1.568, -12.383 },
			coord{ -11.235,   1.210, -10.327 },
			coord{ -11.570,   1.851,  -6.566 },
			prosec_sec_type::BETA_STRAND
		)
		==
		coord{ -11.715646, 0.652680, -9.882344 },
		""
	);
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(prosec_axis_point_gets_alpha_helix_example_correct) {
	static_assert(
		prosec_axis_point_of_residue_triple(
			coord{ -23.119,  12.750, -17.813 },
			coord{ -25.871,  10.535, -16.540 },
			coord{ -23.258,   8.472, -14.671 },
			prosec_sec_type::ALPHA_HELIX
		)
		==
		coord{ -23.585392, 10.599755, -16.286091 },
		""
	);
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
