/// \file
/// \brief The dssp_ss_calc test suite

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

#include "dssp_ss_calc.hpp"

#include <boost/filesystem/path.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/boost_addenda/filesystem/replace_extension_copy.hpp"
#include "common/boost_addenda/test/boost_check_equal_ranges.hpp"
#include "file/dssp_wolf/dssp_file.hpp"
#include "file/dssp_wolf/dssp_file_io.hpp"
#include "file/pdb/pdb.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "structure/sec_struc_calc/dssp/test/dssp_dupl_fixture.hpp"
#include "test/global_test_constants.hpp"

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::sec;
using namespace cath::test;

using boost::filesystem::directory_entry;
using boost::filesystem::directory_iterator;
using boost::filesystem::path;
using boost::range::sort;

namespace cath {
	namespace test {

		struct dssp_ss_calc_test_suite_fixture : protected dssp_dupl_fixture {
		protected:
			~dssp_ss_calc_test_suite_fixture() noexcept = default;

			const path ALIGN_PAIR_SCORE_TEST_DIR{ global_test_constants::TEST_SOURCE_DATA_DIR() / "aligned_pair_score" };

			/// \brief Check the DSSP SS calculations against those in the specified files
			void check_dssp_ss_against_file(const path &arg_dssp_file, ///< The DSSP file to compare against
			                                const path &arg_pdb_file   ///< The PDB file to check
			                                ) {
				const auto parsed_pdb   = read_pdb_file( arg_pdb_file );
				const auto the_protein  = protein_from_dssp_and_pdb( read_dssp_file( arg_dssp_file ), parsed_pdb, dssp_skip_policy::SKIP__BREAK_ANGLES );
				const auto expected_sss = get_sec_strucs( the_protein );

				const auto got_sss      = calc_sec_strucs_of_pdb__recalc_backbone_residues( parsed_pdb );

				BOOST_CHECK_EQUAL_RANGES( got_sss, expected_sss );
			}

			/// \brief Check the DSSP SS calculations against those in the specified file
			///        (assuming the PDB filename is the same as the DSSP's but with the extension stripped off)
			void check_dssp_ss_against_file(const path &arg_dssp_file ///< The DSSP file to compare against
			                                ) {
				check_dssp_ss_against_file( arg_dssp_file, replace_extension_copy( arg_dssp_file ) );
			}

			/// \brief Use the specified directory of DSSP files to check the DSSP SS calculations on the corresponding PDB files
			///        in the specified directory (where the DSSP files are assumed to end ".dssp" and the corresponding PDB files
			///        are assumed to have the same filename but with the extension stripped off)
			void use_dir_of_dssp_files_to_check_hbonds_calcs(const path &arg_dssp_dir, ///< The directory of DSSP files to compare against
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
					// std::cerr << dssp_file << "\n";
					const path pdb_file = replace_extension_copy( arg_pdb_dir / dssp_file.filename() );
					check_dssp_ss_against_file( dssp_file, pdb_file );
				}
			}
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(dssp_ss_calc_test_suite, dssp_ss_calc_test_suite_fixture)


BOOST_AUTO_TEST_CASE(gets_correct_ss_on_1c55A) {
	check_dssp_ss_against_file( ALIGN_PAIR_SCORE_TEST_DIR / "1c55A.dssp" );
}

BOOST_AUTO_TEST_CASE(gets_correct_ss_on_1c56A) {
	check_dssp_ss_against_file( ALIGN_PAIR_SCORE_TEST_DIR / "1c56A.dssp" );
}

BOOST_AUTO_TEST_CASE(gets_correct_ss_on_1hykA) {
	check_dssp_ss_against_file( ALIGN_PAIR_SCORE_TEST_DIR / "1hykA.dssp" );
}

BOOST_AUTO_TEST_CASE(gets_correct_ss_on_1wt7A) {
	check_dssp_ss_against_file( ALIGN_PAIR_SCORE_TEST_DIR / "1wt7A.dssp" );
}

BOOST_AUTO_TEST_CASE(gets_correct_ss_on_example_a) {
	check_dssp_ss_against_file( global_test_constants::EXAMPLE_A_DSSP_FILENAME() );
}

BOOST_AUTO_TEST_CASE(gets_correct_ss_on_example_b) {
	check_dssp_ss_against_file( global_test_constants::EXAMPLE_B_DSSP_FILENAME() );
}



BOOST_AUTO_TEST_SUITE(engineered_examples)



BOOST_AUTO_TEST_CASE(test_without_residues_that_dssp_ignores) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "test_without_residues_that_dssp_ignores.dssp"        );
}

BOOST_AUTO_TEST_CASE(prefer_helix_over_strand) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "prefer_helix_over_strand.dssp"                       );
}



BOOST_AUTO_TEST_SUITE(beta)

BOOST_AUTO_TEST_CASE(beta_bulge_cannot_straddle_break) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "beta_bulge_cannot_straddle_break.dssp"               );
}

BOOST_AUTO_TEST_CASE(beta_bonds_must_exist_from_nh_side) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "beta_bonds_must_exist_from_nh_side.dssp"             );
}

BOOST_AUTO_TEST_CASE(another_beta_bulge) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "another_beta_bulge.dssp"                             );
}

BOOST_AUTO_TEST_CASE(beta_bulges_can_have_equal_gaps) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "beta_bulges_can_have_equal_gaps.dssp"                );
}

BOOST_AUTO_TEST_CASE(ok_consider_beta_bridge_bonds_from_both_sides) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "consider_beta_bridge_bonds_from_both_sides.dssp"     );
}

BOOST_AUTO_TEST_CASE(beta_bridges_forbidden_at_end) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "beta_bridges_forbidden_at_end.dssp"                  );
}

BOOST_AUTO_TEST_CASE(beta_bridge_still_lone_if_neighbours_very_diff) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "beta_bridge_still_lone_if_neighbours_very_diff.dssp" );
}

BOOST_AUTO_TEST_CASE(beta_bonded_residues_must_be_ge_3_apart) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "beta_bonded_residues_must_be_ge_3_apart.dssp"        );
}

BOOST_AUTO_TEST_CASE(do_not_beta_bridge_to_side_of_chain_break) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "do_not_beta_bridge_to_side_of_chain_break.dssp"      );
}

BOOST_AUTO_TEST_CASE(ok_beta_bulge_bonds_can_be_to_same_res) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "beta_bulge_bonds_can_be_to_same_res.dssp"            );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(helix)

BOOST_AUTO_TEST_CASE(helix_bonds_must_exist_from_nh_side) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "helix_bonds_must_exist_from_nh_side.dssp"            );
}

BOOST_AUTO_TEST_CASE(check_helix_bonds_at_both_ends) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "check_helix_bonds_at_both_ends.dssp"                 );
}

BOOST_AUTO_TEST_CASE(not_4_helix_if_5_helix) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "not_4_helix_if_5_helix.dssp"                         );
}

BOOST_AUTO_TEST_CASE(no_5_helix_to_disrupt_4_if_hits_a_3) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "no_5_helix_to_disrupt_4_if_hits_a_3.dssp"            );
}

BOOST_AUTO_TEST_CASE(non_costarting_5_helix_still_disrupts_4_helix) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "non_costarting_5_helix_still_disrupts_4_helix.dssp"  );
}

BOOST_AUTO_TEST_CASE(not_4_if_5_because_no_3_because_4_helix) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "not_4_if_5_because_no_3_because_4_helix.dssp"        );
}

BOOST_AUTO_TEST_CASE(tricky_3_4_5_boundaries) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "tricky_3_4_5_boundaries.dssp"                        );
}

BOOST_AUTO_TEST_CASE(helix_bond_cannot_straddle_break) {
	check_dssp_ss_against_file( DSSP_SS_TEST_DATA_DIR() / "helix_bond_cannot_straddle_break.dssp"               );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE_END()



// BOOST_AUTO_TEST_SUITE(whole_dssp_directory)

// // Testcase for checking all DSSP instances in a directory against the equivalent PDBs in another directory
// BOOST_AUTO_TEST_CASE(whole_dssp_directory_tc) {
// 	use_dir_of_dssp_files_to_check_hbonds_calcs(
// 		"/cath/mothra-data1/people/ucbctnl/dssp-20161211-v4_1_0-after",
// 		"/cath/data/v4_1_0/wholepdb"
// 	);
// }

// BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()
