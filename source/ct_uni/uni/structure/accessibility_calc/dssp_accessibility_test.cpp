/// \file
/// \brief The dssp_accessibility test suite

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

#include "dssp_accessibility.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/size_t_literal.hpp"
#include "file/pdb/pdb.hpp"
#include "structure/geometry/coord.hpp"
#include "test/boost_addenda/boost_check_equal_ranges.hpp"
#include "test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "test/global_test_constants.hpp"

// #include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sec;

using boost::numeric_cast;

BOOST_AUTO_TEST_SUITE(dssp_accessibility_test_suite)

BOOST_AUTO_TEST_CASE(makes_correct_number_of_ball_points) {
	BOOST_CHECK_EQUAL( make_dssp_ball_points( 200 ).size(), 401_z );
}

BOOST_AUTO_TEST_CASE(makes_ball_points_with_unit_y_in_the_middle) {
	BOOST_CHECK_EQUAL( make_dssp_ball_points( 200 )[ 200 ], coord( 0.0, 1.0, 0.0 ) );
}

BOOST_AUTO_TEST_CASE(gets_eg_accessibility_count_correct) {
	const auto parsed_pdb = read_pdb_file( global_test_constants::EXAMPLE_A_PDB_FILENAME() );
	BOOST_REQUIRE_GE   ( parsed_pdb.get_num_residues(), 4_z  );
	const auto &the_res = parsed_pdb.get_residue_of_index__backbone_unchecked( 3 );
	BOOST_REQUIRE_EQUAL( the_res.get_residue_id(), make_residue_id( 'A', 1002 ) );
	BOOST_REQUIRE      ( the_res.has_nitrogen() );
	const pdb_atom &the_atom = the_res.get_nitrogen();
	BOOST_CHECK_EQUAL  ( get_accessibility_count( the_atom, parsed_pdb ), 17_z );
}

BOOST_AUTO_TEST_CASE(gets_eg_accessibility_fraction_correct) {
	const auto parsed_pdb = read_pdb_file( global_test_constants::EXAMPLE_A_PDB_FILENAME() );
	BOOST_REQUIRE_GE   ( parsed_pdb.get_num_residues(), 4_z  );
	const auto &the_res = parsed_pdb.get_residue_of_index__backbone_unchecked( 3 );
	BOOST_REQUIRE_EQUAL( the_res.get_residue_id(), make_residue_id( 'A', 1002 ) );
	BOOST_REQUIRE      ( the_res.has_nitrogen() );
	const pdb_atom &the_atom = the_res.get_nitrogen();
	BOOST_CHECK_EQUAL  ( get_accessibility_fraction( the_atom, parsed_pdb ), 17.0 / 401.0 );
}

BOOST_AUTO_TEST_CASE(gets_eg_accessibility_surface_area_correct) {
	const auto parsed_pdb = read_pdb_file( global_test_constants::EXAMPLE_A_PDB_FILENAME() );
	BOOST_REQUIRE_GE   ( parsed_pdb.get_num_residues(), 4_z  );
	const auto &the_res = parsed_pdb.get_residue_of_index__backbone_unchecked( 3 );
	BOOST_REQUIRE_EQUAL( the_res.get_residue_id(), make_residue_id( 'A', 1002 ) );
	BOOST_REQUIRE      ( the_res.has_nitrogen() );
	const pdb_atom &the_atom = the_res.get_nitrogen();
	BOOST_CHECK_EQUAL  ( get_accessibility_surface_area( the_atom, parsed_pdb ), 4.9558036530705616 );
}

BOOST_AUTO_TEST_CASE(gets_eg_res_accessibility_surface_area_correct) {
	const auto parsed_pdb = read_pdb_file( global_test_constants::EXAMPLE_A_PDB_FILENAME() );
	BOOST_REQUIRE_GE   ( parsed_pdb.get_num_residues(), 4_z  );
	const auto &the_res = parsed_pdb.get_residue_of_index__backbone_unchecked( 3 );
	BOOST_REQUIRE_EQUAL( the_res.get_residue_id(), make_residue_id( 'A', 1002 ) );
	BOOST_CHECK_EQUAL  ( round( get_accessibility_surface_area( the_res, parsed_pdb ) ), 138 );
}

BOOST_AUTO_TEST_CASE(calc_accessibilities_with_scanning_does_not_throw_or_error_on_empty_pdb) {
	BOOST_CHECK_NO_THROW_DIAG( calc_accessibilities_with_scanning( pdb{} ) );
}

BOOST_AUTO_TEST_CASE(gets_pdb_accessibility_surface_area_correct) {
	const auto     parsed_pdb        = read_pdb_file( global_test_constants::EXAMPLE_A_PDB_FILENAME() );
	const size_vec expected_accesses = {
		 231, 150, 190, 138,  58, 177,  74, 137,   2,   3,
		   0,  16,   8,   3,  12,   4,  17,   0,   0,   0,
		   0,   0,   3,  16,   0,  11, 149,  85,  46,  15,
		  50,   4,  45,   0,   1,   5, 120,  14,  10,   9,
		 116,  32,  69,  93,  25, 178,  80,  37,  77,  15,
		  14, 134,  39,  80,  38, 115,  94,  54,  61, 117,
		  87,  79, 127, 131, 118,   0,  81,  84,  13,  22,
		 105,   1,  10, 124,  23,  30, 147, 194, 140,
		   2,  25, 123,  76,   6,  56,  99,   2,  13, 126,
		  40,   0, 131, 145,  81,  51,   9,  76,  59,  71,
		 123,  60,  81,  38,  82,  12,  32,  73,  54,   0,
		   7,  82,  42,   7,  92,   1,   0,   0,  13,  11,
		  40,  41, 169,  39,  37, 190, 106, 291,  86,  43,
		  13,  77,  21,  94,  46,  50,  66,   8,  61,   3,
		  63,  31,  41,  87, 119, 103,  28, 114,  98,  12,
		 100,  68, 200,  31,  84,  48, 202,  82,  68, 210,
		 137,  96, 105,   0,   0,   6,   0,   6,   1,  50,
		  37,  17, 121,  84,   6,  11,  30,   8,   7,   0,
		  52,   1,   0,   0,  84,  42,   0,   0, 107,  52,
		   0,  37,  79,   4,  18, 124, 159, 101,  41, 130,
	};

	// const auto     parsed_pdb        = read_pdb_file( "/cath-tools/refine_stuff/really_long/1u6gC00" );
	// const size_vec expected_accesses = { 1 };

	const size_vec get_accesses_raw  = transform_build<size_vec>(
		calc_accessibilities_with_scanning( parsed_pdb ),
		[] (const double &x) { return numeric_cast<size_t>( round( x ) ); }
	);
	BOOST_CHECK_EQUAL_RANGES( get_accesses_raw, expected_accesses );
}

BOOST_AUTO_TEST_SUITE_END()
