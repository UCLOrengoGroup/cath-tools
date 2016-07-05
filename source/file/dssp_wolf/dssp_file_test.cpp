/// \file
/// \brief The dssp_file test suite

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

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "common/file/simple_file_read_write.h"
#include "common/size_t_literal.h"
#include "file/dssp_wolf/dssp_file.h"
#include "file/dssp_wolf/dssp_file_io.h"
#include "file/dssp_wolf/tally_residue_names.h"
#include "file/dssp_wolf/wolf_file.h"
#include "file/dssp_wolf/wolf_file_io.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_io.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"
#include "test/global_test_constants.h"
#include "test/log_to_ostream_guard.h"

using namespace boost::filesystem;
using namespace boost::math;
using namespace boost::math::constants;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::numeric_cast;
using boost::range::sort;

namespace cath {
	namespace test {

		/// \brief The dssp_wolf_file_test_suite_fixture to assist in testing dssp_file
		struct dssp_wolf_file_test_suite_fixture : protected global_test_constants {
		protected:
			~dssp_wolf_file_test_suite_fixture() noexcept = default;

			const path wolf_file_of_example_id(const string &);
			const path dssp_file_of_example_id(const string &);
			const path pdb_file_of_example_id(const string &);
			const path expected_wolf_phi_psi_of_example_id(const string &);

			void check_phi_psi_from_wolf_against_expected(const string &);
			void compare_phi_psi_from_pdb_vs_dssp(const string &);
			void check_pdb_and_dssp_built_protein(const string &);
			void compare_residue_frames_from_pdb_vs_wolf(const string &);

			const string EXAMPLE_WOLF_FILE_SUFFIX     { ".wolf"                  };
			const string EXAMPLE_DSSP_FILE_SUFFIX     { ".dssp"                  };
			const string EXAMPLE_PDB_FILE_SUFFIX      { ""                       };
			const string EXPECTED_WOLF_PHI_PSI_SUFFIX { ".expected_wolf_phi_psi" };
		};

	}
}

/// \brief TODOCUMENT
const path cath::test::dssp_wolf_file_test_suite_fixture::wolf_file_of_example_id(const string &arg_example_id ///< TODOCUMENT
                                                                                  ) {
	return TEST_RESIDUE_NAMES_DATA_DIR() / (arg_example_id + EXAMPLE_WOLF_FILE_SUFFIX);
}

/// \brief TODOCUMENT
const path cath::test::dssp_wolf_file_test_suite_fixture::dssp_file_of_example_id(const string &arg_example_id ///< TODOCUMENT
                                                                                  ) {
	return TEST_RESIDUE_NAMES_DATA_DIR() / (arg_example_id + EXAMPLE_DSSP_FILE_SUFFIX);
}

/// \brief TODOCUMENT
const path cath::test::dssp_wolf_file_test_suite_fixture::pdb_file_of_example_id(const string &arg_example_id ///< TODOCUMENT
                                                                                 ) {
	return TEST_RESIDUE_NAMES_DATA_DIR() / (arg_example_id + EXAMPLE_PDB_FILE_SUFFIX);
}

/// \brief TODOCUMENT
const path cath::test::dssp_wolf_file_test_suite_fixture::expected_wolf_phi_psi_of_example_id(const string &arg_example_id ///< TODOCUMENT
                                                                                              ) {
	return TEST_RESIDUE_NAMES_DATA_DIR() / (arg_example_id + EXPECTED_WOLF_PHI_PSI_SUFFIX);
}

/// \brief TODOCUMENT
void cath::test::dssp_wolf_file_test_suite_fixture::check_phi_psi_from_wolf_against_expected(const string &arg_example_id ///< TODOCUMENT
                                                                                             ) {
	// Get the phi and psi angles from the wolf file
	const protein     wolf_protein = protein_from_wolf( read_wolf( wolf_file_of_example_id(arg_example_id) ) );
	doub_vec got_wolf_phi_psi_values;
	for (const residue &the_residue : wolf_protein) {
		got_wolf_phi_psi_values.push_back( round( angle_in_degrees( the_residue.get_phi_angle() ) ) );
		got_wolf_phi_psi_values.push_back( round( angle_in_degrees( the_residue.get_psi_angle() ) ) );
	}

//	// To update file, uncomment this, run all tests, then execute:
//	//    mv build-test-data/residue_names/id.new_got_wolf_phi_psi build-test-data/residue_names/id.expected_wolf_phi_psi
//	write_file(
//		TEST_RESIDUE_NAMES_DATA_DIR() / ( arg_example_id + ".new_got_wolf_phi_psi"),
//		got_wolf_phi_psi_values
//	);

	// Load the expected phi and psi angles and compare them to what we actually got
	const doub_vec expected_wolf_phi_psi_values = read_file<double>( expected_wolf_phi_psi_of_example_id( arg_example_id ) );
	BOOST_CHECK_EQUAL_RANGES( expected_wolf_phi_psi_values, got_wolf_phi_psi_values );
}

/// \brief TODOCUMENT
///
/// The NULL records in the DSSP files make this a bit trickier.
/// Should ideally do an NW over the residue names here but the
/// NW code isn't in shape for that use yet so this is a bit hackier.
void cath::test::dssp_wolf_file_test_suite_fixture::compare_phi_psi_from_pdb_vs_dssp(const string &arg_example_id ///< TODOCUMENT
                                                                                     ) {
	const auto the_pdb_file      = read_pdb_file(   pdb_file_of_example_id( arg_example_id ) );
	const auto the_dssp_file     = read_dssp_file( dssp_file_of_example_id( arg_example_id ) );
	const auto pdb_prot          = build_protein_of_pdb( the_pdb_file      );
	const auto num_dssp_residues = the_dssp_file.get_num_residues();
	const auto num_pdb_residues  = pdb_prot.get_length();

	// Check that the number of residues in each matches to within 10%
	BOOST_REQUIRE_CLOSE(
		numeric_cast<double>( num_dssp_residues ),
		numeric_cast<double>( num_pdb_residues  ),
		10.0
	);

	// Grab the residues names from the DSSP and PDB and then tally them up
	const auto pdb_res_names  = the_pdb_file.get_residue_names_of_first_chain__backbone_unchecked();
	const auto dssp_res_names = get_residue_names(the_dssp_file, false);

	const auto alignment      = tally_residue_names(pdb_res_names, dssp_res_names, false);

	// Loop over the tallied residue names
	for (const size_size_pair &aligned_pair : alignment) {
		// Grab the indices
		const size_t  &pdb_res_index    = aligned_pair.first;
		const size_t  &dssp_res_index   = aligned_pair.second;

		// Grab the residues
		const residue &dssp_residue     = the_dssp_file.get_residue_of_index(dssp_res_index);
		const residue &pdb_prot_residue = pdb_prot.get_residue_ref_of_index( pdb_res_index );

		BOOST_REQUIRE_EQUAL( pdb_prot_residue.get_pdb_residue_name(), dssp_residue.get_pdb_residue_name() );

		// The DSSP phi/psi angles get rounded to one decimal place and then an integer
		// which means that, for example, 10.47 gets rounded to 11
		// so the test needs to do the same to match
		const doub_angle dssp_phi = dssp_residue.get_phi_angle();
		const doub_angle dssp_psi = dssp_residue.get_psi_angle();
		if ( dssp_phi != residue::DEFAULT_PHI_PSI() ) {
			BOOST_CHECK_CLOSE( angle_in_degrees( dssp_phi ), round( round( 10.0 * angle_in_degrees( pdb_prot_residue.get_phi_angle() ) ) / 10.0 ), ACCURACY_PERCENTAGE() );
		}
		if ( dssp_psi != residue::DEFAULT_PHI_PSI() ) {
			BOOST_CHECK_CLOSE( angle_in_degrees( dssp_psi ), round( round( 10.0 * angle_in_degrees( pdb_prot_residue.get_psi_angle() ) ) / 10.0 ), ACCURACY_PERCENTAGE() );
		}
	}
}

/// \brief TODOCUMENT
void cath::test::dssp_wolf_file_test_suite_fixture::check_pdb_and_dssp_built_protein(const string &arg_example_id ///< TODOCUMENT
                                                                                     ) {
	const auto the_pdb_file               = read_pdb_file (  pdb_file_of_example_id( arg_example_id ) );
	const auto the_dssp_file              = read_dssp_file( dssp_file_of_example_id( arg_example_id ) );
	const auto pdb_prot                   = build_protein_of_pdb( the_pdb_file      );
	const auto num_non_null_dssp_residues = get_num_non_null_residues( the_dssp_file );
	const auto num_pdb_residues           = pdb_prot.get_length();

	ostringstream test_ss;
	const log_to_ostream_guard the_guard{ test_ss };

	const protein      combi_prot_with_all_pdb_residues   = protein_from_dssp_and_pdb(the_dssp_file, the_pdb_file, false );
	const protein      combi_prot_with_dssp_only_residues = protein_from_dssp_and_pdb(the_dssp_file, the_pdb_file, true  );

	BOOST_REQUIRE_EQUAL( num_non_null_dssp_residues, combi_prot_with_dssp_only_residues.get_length() );
	BOOST_REQUIRE_EQUAL( num_pdb_residues,           combi_prot_with_all_pdb_residues.get_length()   );

	// Compare the combi's residues with the DSSP residues
	//( this requires two counters to handle skipping null DSSP residues)
	size_t combi_residue_ctr = 0;
	for (const residue dssp_residue : the_dssp_file) {

		// If this is a null DSSP residue then just move onto the next one
		if ( is_null_residue(dssp_residue) ) {
			continue;
		}
		const residue built_residue = combi_prot_with_dssp_only_residues.get_residue_ref_of_index( combi_residue_ctr );
		BOOST_CHECK_EQUAL( dssp_residue.get_pdb_residue_name(),   built_residue.get_pdb_residue_name()   );
		BOOST_CHECK_EQUAL( dssp_residue.get_sec_struc_number(),   built_residue.get_sec_struc_number()   );
		BOOST_CHECK_EQUAL( dssp_residue.get_sec_struc_type(),     built_residue.get_sec_struc_type()     );
//		cerr << dssp_residue  << endl;
//		cerr << built_residue << endl;
//		cerr << endl;
		++combi_residue_ctr;
	}

	for (size_t pdb_residue_ctr = 0; pdb_residue_ctr < num_pdb_residues; ++pdb_residue_ctr) {
		const pdb_residue the_pdb_residue = the_pdb_file.get_residue_cref_of_index__backbone_unchecked(pdb_residue_ctr);
		const residue     built_residue   = combi_prot_with_all_pdb_residues.get_residue_ref_of_index( pdb_residue_ctr );
		BOOST_CHECK_EQUAL( the_pdb_residue.get_residue_name(),                             built_residue.get_pdb_residue_name()  );
		BOOST_CHECK_EQUAL( get_carbon_alpha_coord_of_residue( the_pdb_residue ),           built_residue.get_carbon_alpha_coord() );
		BOOST_CHECK_EQUAL( get_or_predict_carbon_beta_coord_of_residue( the_pdb_residue ), built_residue.get_carbon_beta_coord()  );
//		cerr << pdb_residue   << endl;
//		cerr << built_residue << endl;
//		cerr << endl;
	}
}

/// \brief TODOCUMENT
void cath::test::dssp_wolf_file_test_suite_fixture::compare_residue_frames_from_pdb_vs_wolf(const string &arg_example_id ///< TODOCUMENT
                                                                                            ) {
	const wolf_file            the_wolf_file     = read_wolf(     wolf_file_of_example_id( arg_example_id ) );
	const pdb                  the_pdb_file      = read_pdb_file(  pdb_file_of_example_id( arg_example_id ) );
	const protein              pdb_prot          = build_protein_of_pdb( the_pdb_file      );
	const wolf_file::size_type num_wolf_residues = the_wolf_file.get_num_residues();
	const size_t   num_pdb_residues  = pdb_prot.get_length();

	BOOST_REQUIRE_EQUAL(num_pdb_residues, num_wolf_residues);

	doub_angle_vec angles_between_frame_pairs;
	for (size_t residue_ctr = 0; residue_ctr < num_pdb_residues; ++residue_ctr) {
		const residue &wolf_residue     = the_wolf_file.get_residue_of_index( residue_ctr );
		const residue &pdb_prot_residue = pdb_prot.get_residue_ref_of_index( residue_ctr );

		BOOST_REQUIRE_EQUAL( pdb_prot_residue.get_pdb_residue_name(), wolf_residue.get_pdb_residue_name() );

		const rotation   raw_wolf_frame   = wolf_residue.get_frame();
		const rotation   wolf_frame       = tidy_copy(raw_wolf_frame, 0.001);
		const rotation   pdb_frame        = pdb_prot_residue.get_frame();
		const doub_angle the_angle        = angle_between_rotations( wolf_frame, pdb_frame );
		angles_between_frame_pairs.push_back( the_angle );
	}

	BOOST_REQUIRE_GT( angles_between_frame_pairs.size(), 0_z );

	// Sort the list of angles in ascending order
	sort( angles_between_frame_pairs );

	// Check that around 95% of the way through the list of angles, the values are still
	// less than 3 degrees (in radians)
	const size_t ninety_fifth_perc_index = numeric_cast<size_t>(
		round( 0.95 * numeric_cast<double>( angles_between_frame_pairs.size() ) )
	);
	BOOST_CHECK_LT(
		angles_between_frame_pairs[ ninety_fifth_perc_index ],
		make_angle_from_degrees<double>( 3.0 )
	);
}
/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(dssp_wolf_file_test_suite, cath::test::dssp_wolf_file_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(compare_phi_psi_from_pdb_vs_dssp__1cukA) {
	compare_phi_psi_from_pdb_vs_dssp("1cukA");
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(compare_phi_psi_from_pdb_vs_dssp__1bmv2) {
	compare_phi_psi_from_pdb_vs_dssp("1bmv2");
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(compare_phi_psi_from_pdb_vs_dssp__1p3cA) {
	compare_phi_psi_from_pdb_vs_dssp("1p3cA");
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(check_pdb_and_dssp_built_protein__1cukA) {
	check_pdb_and_dssp_built_protein("1cukA");
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(check_pdb_and_dssp_built_protein__1bmv2) {
	check_pdb_and_dssp_built_protein("1bmv2");
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(check_pdb_and_dssp_built_protein__1p3cA) {
	check_pdb_and_dssp_built_protein("1p3cA");
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(compare_residue_frames_from_pdb_vs_wolf__1bmv2) {
	compare_residue_frames_from_pdb_vs_wolf("1bmv2");
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(insert_code) {
	const size_chain_residue_tuple parsed_line = parse_dssp_residue_line(
		"  250  918AA T              0   0  200     -2,-0.2    -1,-0.2    -3,-0.0    -2,-0.1   0.692 360.0 360.0-115.3 360.0   -7.3   40.4   19.8"
	);
	BOOST_CHECK_EQUAL(                   250_z,  get<0>(parsed_line )                        );
	BOOST_CHECK_EQUAL(      chain_label( 'A' ),  get<1>(parsed_line )                        );
	BOOST_CHECK_EQUAL( residue_name(918, 'A') ,  get<2>(parsed_line ).get_pdb_residue_name() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(null_dssp_residue) {
	const size_chain_residue_tuple parsed_line = parse_dssp_residue_line(
		"   36        !              0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0"
	);
	BOOST_CHECK_EQUAL(                    0_z, get<0>( parsed_line ) );
	BOOST_CHECK_EQUAL(    chain_label( ' ' ), get<1>( parsed_line ) );
	BOOST_CHECK_EQUAL( residue::NULL_RESIDUE, get<2>( parsed_line ) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(wolf_phi_psi__1a04A02) {
	check_phi_psi_from_wolf_against_expected( "1a04A02" );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(wolf_phi_psi__1fseB00) {
	check_phi_psi_from_wolf_against_expected( "1fseB00" );
}

BOOST_AUTO_TEST_SUITE_END()
