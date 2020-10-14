/// \file
/// \brief The sec_calc test suite

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

#include "sec_calc.hpp"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/test/unit_test.hpp>

#include "common/boost_addenda/range/front.hpp"
#include "common/pair_insertion_operator.hpp"
#include "file/dssp_wolf/dssp_file.hpp"
#include "file/dssp_wolf/dssp_file_io.hpp"
#include "file/pdb/pdb.hpp"
#include "file/sec/sec_file.hpp"
#include "file/sec/sec_file_io.hpp"
#include "file/sec/sec_file_record.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "test/boost_addenda/boost_check_equal_ranges.hpp"
#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sec;

BOOST_FIXTURE_TEST_SUITE(sec_calc_test_suite, global_test_constants)

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

BOOST_AUTO_TEST_CASE(get_sec_starts_and_stops_works) {
	const protein the_protein = protein_from_dssp_and_pdb(
		read_dssp_file( EXAMPLE_B_DSSP_FILENAME() ),
		read_pdb_file ( EXAMPLE_B_PDB_FILENAME()  )
	);
	const auto got_starts_stops = get_sec_starts_and_stops( the_protein );
	const size_size_pair_vec expected_starts_stops = {
		{   4,   8 },
		{  13,  24 },
		{  28,  33 },
		{  48,  51 },
		{  57,  64 },
		{  69,  72 },
		{  85, 100 },
		{ 104, 108 },
		{ 125, 140 },
		{ 144, 148 },
		// { 152, 154 }, // The sec file differs here, it has: { 151, 155 },
		// { 163, 165 }, // The sec file differs here, it has: { 162, 166 },
		// { 174, 176 }, // The sec file differs here, it has: { 173, 177 },
		{ 177, 186 },
		{ 197, 201 },
	};
	BOOST_CHECK_EQUAL_RANGES( got_starts_stops, expected_starts_stops );
}

BOOST_AUTO_TEST_CASE(pymol_of_get_sec_records_works) {
	const protein the_protein = protein_from_dssp_and_pdb(
		read_dssp_file( EXAMPLE_B_DSSP_FILENAME() ),
		read_pdb_file ( EXAMPLE_B_PDB_FILENAME()  )
	);
	BOOST_CHECK( boost::algorithm::contains( get_pymol_script_text( the_protein, get_sec_records( the_protein ) ), "show_as" ) );
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
