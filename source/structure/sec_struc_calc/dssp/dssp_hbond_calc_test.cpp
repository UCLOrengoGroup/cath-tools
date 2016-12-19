/// \file
/// \brief The dssp_hbond_calc test suite

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

#include "dssp_hbond_calc.hpp"

#include <boost/filesystem.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "common/boost_addenda/filesystem/replace_extension_copy.hpp"
#include "common/size_t_literal.hpp"
#include "file/dssp_wolf/dssp_file.hpp"
#include "file/dssp_wolf/dssp_file_io.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "structure/sec_struc_calc/dssp/bifur_hbond_list.hpp"
#include "test/dssp/dssp_dupl_fixture.hpp"
#include "test/global_test_constants.hpp"
// #include "common/chrono/duration_to_seconds_string.hpp"

// #include <chrono>

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sec;
using namespace cath::test;

using boost::filesystem::directory_entry;
using boost::filesystem::directory_iterator;
using boost::filesystem::path;
using boost::irange;
using boost::none;
using boost::range::sort;
// using std::chrono::high_resolution_clock;

namespace cath {
	namespace test {

		/// \brief The dssp_hbond_calc_test_suite_fixture to assist in testing dssp_hbond_calc
		struct dssp_hbond_calc_test_suite_fixture : protected dssp_dupl_fixture {
		protected:
			~dssp_hbond_calc_test_suite_fixture() noexcept = default;

			void use_dssp_file_to_check_hbonds_calcs(const path &,
			                                         const path &,
			                                         const bool & = false);

			void use_dssp_file_to_check_hbonds_calcs(const path &,
			                                         const bool & = false);

			void use_dir_of_dssp_files_to_check_hbonds_calcs(const path &,
			                                                 const path &);
		};

		/// \brief Use the specified DSSP file to check the hbond calculations for the specified PDB file
		///
		/// \todo Remove the arg_warn_only_on_diff argument if dssp_ignores_valid_residue_1 is resolved
		inline void dssp_hbond_calc_test_suite_fixture::use_dssp_file_to_check_hbonds_calcs(const path &arg_dssp_file,        ///< The DSSP file to compare against
		                                                                                    const path &arg_pdb_file,         ///< The PDB file to check
		                                                                                    const bool &arg_warn_only_on_diff ///< Whether to only warn on differences
		                                                                                    ) {
			BOOST_REQUIRE( exists( arg_pdb_file ) );
			// const auto read_dssp_start_time = high_resolution_clock::now();
			const auto dssp_hbonds = parse_dssp_for_calc_testing( arg_dssp_file );
			// const auto read_dssp_stop_time  = high_resolution_clock::now();

			if ( ! dssp_hbonds.empty() ) {
				// const auto read_pdb_start_time = high_resolution_clock::now();
				const auto parsed_pdb          = read_pdb_file( arg_pdb_file );
				// const auto read_pdb_stop_time  = high_resolution_clock::now();

				protein_from_dssp_and_pdb( read_dssp_file( arg_dssp_file ), parsed_pdb, dssp_skip_policy::DONT_SKIP__BREAK_ANGLES );

				// const auto calc_start_time = high_resolution_clock::now();
				const auto bifur_hbonds = dssp_hbond_calc::calc_bifur_hbonds_of_pdb__recalc_backbone_residues(
					parsed_pdb
				);
				// const auto calc_stop_time  = high_resolution_clock::now();

				// std::cerr << "Read DSSP   : " << durn_to_seconds_string ( read_dssp_stop_time - read_dssp_start_time ) << "\n";
				// std::cerr << "Read PDB    : " << durn_to_seconds_string ( read_pdb_stop_time  - read_pdb_start_time  ) << "\n";
				// std::cerr << "Calc hbonds : " << durn_to_seconds_string ( calc_stop_time      - calc_start_time      ) << "\n";

				// BOOST_TEST_INFO() isn't present in Boost > 1.58.0
				// BOOST_TEST_INFO  ( "Checking DSSP file \"" + arg_dssp_file.string() + "\"" );
				if ( arg_warn_only_on_diff ) {
					BOOST_WARN_EQUAL ( difference_string( dssp_hbonds, bifur_hbonds ), none );
				}
				else {
					BOOST_CHECK_EQUAL( difference_string( dssp_hbonds, bifur_hbonds ), none );
				}
			}
		}

		/// \brief Use the specified DSSP file to check the hbond calculations for the PDB file with the same filename but with the extension stripped off
		///
		/// \todo Remove the arg_warn_only_on_diff argument if dssp_ignores_valid_residue_1 is resolved
		inline void dssp_hbond_calc_test_suite_fixture::use_dssp_file_to_check_hbonds_calcs(const path &arg_dssp_file,        ///< The DSSP file to compare against
		                                                                                    const bool &arg_warn_only_on_diff ///< Whether to only warn on differences
		                                                                                    ) {
			use_dssp_file_to_check_hbonds_calcs( arg_dssp_file, replace_extension_copy( arg_dssp_file ), arg_warn_only_on_diff );
		}

		/// \brief Use the specified directory of DSSP files to check the hbond calculations on the corresponding PDB files
		///        in the specified directory (where the DSSP files are assumed to end ".dssp" and the corresponding PDB files
		///        are assumed to have the same filename but with the extension stripped off)
		inline void dssp_hbond_calc_test_suite_fixture::use_dir_of_dssp_files_to_check_hbonds_calcs(const path &arg_dssp_dir, ///< The directory of DSSP files to compare against
		                                                                                            const path &arg_pdb_dir   ///< The directory of PDB files to check
		                                                                                            ) {
			BOOST_REQUIRE( is_directory( arg_dssp_dir ) );
			BOOST_REQUIRE( is_directory( arg_pdb_dir  ) );

			const path_vec sorted_dssp_files = [&] () {
				path_vec results;
				for (const directory_entry &dssp_entry : directory_iterator( arg_dssp_dir ) ) {
					const path &dssp_file = dssp_entry.path();
					if ( boost::ends_with( dssp_file.string(), ".dssp" ) ) {
						results.push_back( dssp_file );
					}
				}
				sort( results );
				return results;
			} ();

			for (const path &dssp_file : sorted_dssp_files ) {
				const path pdb_file = replace_extension_copy( arg_pdb_dir / dssp_file.filename() );
				use_dssp_file_to_check_hbonds_calcs( dssp_file, pdb_file );
			}
		}

	}
}

BOOST_FIXTURE_TEST_SUITE(dssp_hbond_calc_test_suite, dssp_hbond_calc_test_suite_fixture)


BOOST_AUTO_TEST_SUITE(issues)

/// \brief This appears to be a minor error in DSSP that can cause a valid residue (eg A124 in 1sx2) to be ignored
///        depending on the residue that precedes it. I've contacted the DSSP people and we hope to resolve this in DSSP.
BOOST_AUTO_TEST_CASE(dssp_ignores_valid_residue_1) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "dssp_ignores_valid_residue_1.dssp", true );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(engineered_test_examples)


BOOST_AUTO_TEST_CASE(heeds_ter_record_blocking_bonds) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "ter_record_blocks_bonds_1.dssp"   );
}

BOOST_AUTO_TEST_CASE(rejects_residue_with_nonstd_first_altlocn) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "nonstd_first_altlocn_1.dssp"      );
}

BOOST_AUTO_TEST_CASE(prevents_nhbond_for_first_residue_after_chain_break) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "chain_break_1.dssp"               );
}

BOOST_AUTO_TEST_CASE(sets_energy_to_minimum_when_very_close) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "very_close_1.dssp"                );
}

BOOST_AUTO_TEST_CASE(takes_last_of_each_atom_for_multifold_residues) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "multifold_residue_1.dssp"         );
}

BOOST_AUTO_TEST_CASE(works_with_non_std_amino_acid_hetatm) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "non_std_amino_acid_hetatm_1.dssp" );
}

BOOST_AUTO_TEST_CASE(rounds_energy_value_like_dssp_1) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "rounding_1.dssp"                  );
}

BOOST_AUTO_TEST_CASE(rounds_energy_value_like_dssp_2) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "rounding_2.dssp"                  );
}

BOOST_AUTO_TEST_CASE(prefers_earlier_bonds_1) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "prefer_earlier_bonds_1.dssp"      );
}

BOOST_AUTO_TEST_CASE(calculates_distance_within_cutoff_using_double_precision) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "bonded_residues_almost_nine_angstroms_apart.dssp" );
}

BOOST_AUTO_TEST_CASE(prob_hbond_to_residue_at_start_of_chain_1) {
	use_dssp_file_to_check_hbonds_calcs( DSSP_HBOND_TEST_DATA_DIR() / "hbond_to_residue_at_start_of_chain_1.dssp"        );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(complete_dssp_files)

BOOST_AUTO_TEST_CASE(duplicates_dssp_hbonds_on_normal) {
	use_dssp_file_to_check_hbonds_calcs( global_test_constants::EXAMPLE_A_DSSP_FILENAME() );
}

BOOST_AUTO_TEST_CASE(handles_mixed_alt_loc) {
	use_dssp_file_to_check_hbonds_calcs( global_test_constants::EXAMPLE_B_DSSP_FILENAME() );
}

BOOST_AUTO_TEST_CASE(handles_residues_dssp_skips_due_to_alt_locs) {
	use_dssp_file_to_check_hbonds_calcs( global_test_constants::TEST_RESIDUE_IDS_DATA_DIR() / "1p3cA.dssp" ); // Eg, residue 36
}

BOOST_AUTO_TEST_CASE(handles_incomplete_residues) {
	use_dssp_file_to_check_hbonds_calcs( global_test_constants::TEST_SOURCE_DATA_DIR() / "incomplete_residues" / "1y1vS.dssp" );
}

BOOST_AUTO_TEST_SUITE_END()



// BOOST_AUTO_TEST_SUITE(whole_dssp_directory)

// // Testcase for checking all DSSP instances in a directory against the equivalent PDBs in another directory
// BOOST_AUTO_TEST_CASE(whole_dssp_directory_tc) {
// 	use_dir_of_dssp_files_to_check_hbonds_calcs(
// 		"/cath/mothra-data1/people/ucbctnl/temp_tony_20160601/dssp",
// 		"/cath/mothra-data1/people/ucbctnl/temp_tony_20160601/pdb"

// 		// "/cath/data/v4_1_0/wholedssp/",
// 		// "/cath/data/v4_1_0/wholepdb/"

// 		// "/cath-tools/build-test-data/residue_names",
// 		// "/cath-tools/build-test-data/residue_names"
// 	);
// }

// BOOST_AUTO_TEST_SUITE_END()



// BOOST_AUTO_TEST_SUITE(speed_test)

// BOOST_AUTO_TEST_CASE(speed_test_tc) {
// 	use_dssp_file_to_check_hbonds_calcs( "/cath-tools/dssp_stuff/1mvw.dssp" );
// }

// BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE_END()
