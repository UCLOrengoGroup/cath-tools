/// \file
/// \brief The ssap test suite

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

#include <boost/optional.hpp> // ***** TEMPORARY *****

#include "chopping/domain/domain.hpp"
#include "chopping/region/region.hpp"
#include "common/boost_addenda/log/log_to_ostream_guard.hpp"
#include "common/boost_addenda/test/boost_check_equal_ranges.hpp"
#include "common/file/simple_file_read_write.hpp"
#include "common/size_t_literal.hpp"
#include "common/type_aliases.hpp"
#include "file/options/data_dirs_options_block.hpp"
#include "ssap/ssap.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_source_file_set/protein_from_wolf_and_sec.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::none;

namespace cath {
	namespace test {

		/// \brief The ssap_test_suite_fixture to assist in testing ssap functions
		///
		/// Note: it's important that this fixture inherits from global_test_constants so that the
		/// setup for these tests calls reset_ssap_global_variables()
		struct ssap_test_suite_fixture : protected global_test_constants {
		protected:
			ssap_test_suite_fixture() {
				reset_ssap_global_variables();
			}
			~ssap_test_suite_fixture() noexcept = default;

		public:
			static const string id_1a04A02;
			static const string id_1fseB00;
		};

	}  // namespace test
}  // namespace cath

/// \brief A string containing the ID for 1a04A02 to be used in a ssap_pair_fixture
const string cath::test::ssap_test_suite_fixture::id_1a04A02("1a04A02");

/// \brief A string containing the ID for 1fseB00 to be used in a ssap_pair_fixture
const string cath::test::ssap_test_suite_fixture::id_1fseB00("1fseB00");

namespace cath {
	namespace test {

		/// \brief A fixture for testing specific SSAP pairs
		///
		/// At present, this just parses the files for the pair on setup and
		/// creates a const-reference alias for each resulting protein
		template <
			const string * const ID1,
			const string * const ID2
		>
		struct ssap_pair_fixture : protected ssap_test_suite_fixture {
		protected:
			~ssap_pair_fixture() noexcept = default;

		public:
			ostringstream log_output_ss;
			log_to_ostream_guard log_output_guard{ log_output_ss };
			const string id1 = { *ID1 };
			const string id2 = { *ID2 };
			const data_dirs_spec data_dirs       = build_data_dirs_spec_of_dir( TEST_SSAP_REGRESSION_DATA_DIR() );
			const prot_prot_pair parsed_proteins = read_protein_pair( id1, none, id2, none, data_dirs, protein_from_wolf_and_sec(), none );
			const protein &prot1 = parsed_proteins.first;
			const protein &prot2 = parsed_proteins.second;
			
			void check_context_sec_scores_as_expected() const;

			void check_residues_have_similar_area_angle_props() const;
		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
template < const string * const ID1, const string * const ID2 >
void cath::test::ssap_pair_fixture<ID1, ID2>::check_context_sec_scores_as_expected() const {
	const size_t num_sec_strucs_1 = prot1.get_num_sec_strucs();
	const size_t num_sec_strucs_2 = prot2.get_num_sec_strucs();
	vector<score_type> got_context_sec_scores;
	got_context_sec_scores.reserve(
		  num_sec_strucs_1 * num_sec_strucs_1
		* num_sec_strucs_2 * num_sec_strucs_2
	);
	for (size_t sec_struc_from_ctr_1 = 0; sec_struc_from_ctr_1 < num_sec_strucs_1; ++sec_struc_from_ctr_1) {
		for (size_t sec_struc_from_ctr_2 = 0; sec_struc_from_ctr_2 < num_sec_strucs_2; ++sec_struc_from_ctr_2) {
			for (size_t sec_struc_to_ctr_1 = 0; sec_struc_to_ctr_1 < num_sec_strucs_1; ++sec_struc_to_ctr_1) {
				for (size_t sec_struc_to_ctr_2 = 0; sec_struc_to_ctr_2 < num_sec_strucs_2; ++sec_struc_to_ctr_2) {
					const score_type context_sec_score = context_sec(
						prot1,
						prot2,
						sec_struc_from_ctr_1,
						sec_struc_from_ctr_2,
						sec_struc_to_ctr_1,
						sec_struc_to_ctr_2
					);
					got_context_sec_scores.push_back(context_sec_score);
				}
			}
		}
	}

//	cerr << endl << endl << "SCORES:" << endl;
//	for (const context_sec_type &score : got_context_sec_scores) {
//		cerr << " " << score;
//	}
//	cerr << endl << endl;
	const vector<score_type> expected_context_sec_scores = read_file<score_type>(
		TEST_SSAP_REGRESSION_DATA_DIR() / ( id1 + "_" + id2 + ".expected_context_sec_scores" )
	);
	BOOST_CHECK_EQUAL_RANGES( expected_context_sec_scores, got_context_sec_scores );
}

/// \brief TODOCUMENT
template < const string * const ID1, const string * const ID2 >
void cath::test::ssap_pair_fixture<ID1, ID2>::check_residues_have_similar_area_angle_props() const {
	const size_t num_residues_1 = prot1.get_length();
	const size_t num_residues_2 = prot2.get_length();
	vector<bool> got_residues_similar;
	got_residues_similar.reserve( num_residues_1 * num_residues_2 );
	for (size_t residue_ctr_1 = 0; residue_ctr_1 < num_residues_1; ++residue_ctr_1) {
		for (size_t residue_ctr_2 = 0; residue_ctr_2 < num_residues_2; ++residue_ctr_2) {
			const bool residues_similar = residues_have_similar_area_angle_props(
				prot1.get_residue_ref_of_index( residue_ctr_1 ),
				prot2.get_residue_ref_of_index( residue_ctr_2 )
			);
			got_residues_similar.push_back( residues_similar );
		}
	}

//	// To update file, uncomment this, run all tests, then execute:
//	//    mv build-test-data/ssap_regression/id1_id2.new_got_residues_similar build-test-data/ssap_regression/id1_id2.expected_residues_similar
//	write_file(
//		TEST_SSAP_REGRESSION_DATA_DIR() / ( id1 + "_" + id2 + ".new_got_residues_similar"),
//		got_residues_similar
//	);

	const vector<bool> expected_residues_similar = read_file<bool>(
		TEST_SSAP_REGRESSION_DATA_DIR() / ( id1 + "_" + id2 + ".expected_residues_similar" )
	);
	BOOST_CHECK_EQUAL_RANGES( expected_residues_similar, got_residues_similar );
}

/// \todo Should add further regression tests (not least for context_res() )
//
//int context_res(const residue &,
//                const residue &,
//                const residue &,
//                const residue &);

BOOST_FIXTURE_TEST_SUITE(ssap_test_suite, cath::test::ssap_test_suite_fixture)

/// \brief Type alias for a fixture for testing the pair 1a04A02 / 1fseB00
///
/// 1a04A02 / 1fseB00 is a good pair for regression testing SSAPs for a few reasons:
///  * both fairly small (80/70 residues respectively)
///  * both have multiple secondary structures in sec files (5/4 respectively)
///  * good SSAP score (87.49) so the answer shouldn't be very ambiguous
///  * they exhibited the regression being investigated at the time of writing
using fixture_1a04A02_1fseB00 = cath::test::ssap_pair_fixture<&cath::test::ssap_test_suite_fixture::id_1a04A02,
                                                              &cath::test::ssap_test_suite_fixture::id_1fseB00> ;

/// \brief Pair of tests to check that both can expect an example SSAP global variable (global_run_counter)
///        to be reset, even if both tests alter it before finishing.
///
/// The two tests perform both roles so that it doesn't matter if their order is randomised
BOOST_FIXTURE_TEST_CASE(check_globals_are_being_reset_before_tests_1, fixture_1a04A02_1fseB00) {
	BOOST_CHECK_EQUAL(0, temp_get_global_run_counter());
	temp_set_global_run_counter(1234);
}

/// \brief Pair of tests to check that both can expect an example SSAP global variable (global_run_counter)
///        to be reset, even if both tests alter it before finishing.
///
/// The two tests perform both roles so that it doesn't matter if their order is randomised
BOOST_FIXTURE_TEST_CASE(check_globals_are_being_reset_before_tests_2, fixture_1a04A02_1fseB00) {
	BOOST_CHECK_EQUAL(0, temp_get_global_run_counter());
	temp_set_global_run_counter(1234);
}

/// \brief Check that 1a04A02 has 5 secondary structures
BOOST_FIXTURE_TEST_CASE(prot_1a04A02_has_5_sec_strucs, fixture_1a04A02_1fseB00) {
	BOOST_CHECK_EQUAL(prot1.get_num_sec_strucs(), 5_z); // 1a04A02
}

/// \brief Check that 1fseB00 has 4 secondary structures
BOOST_FIXTURE_TEST_CASE(prot_1fseB00_has_4_sec_strucs, fixture_1a04A02_1fseB00) {
	BOOST_CHECK_EQUAL(prot2.get_num_sec_strucs(), 4_z); // 1fseB00
}

/// \brief Check that 1a04A02 has 80 residues
BOOST_FIXTURE_TEST_CASE(prot_1a04A02_has_80_residues, fixture_1a04A02_1fseB00) {
	BOOST_CHECK_EQUAL( prot1.get_length(), 80_z ); // 1a04A02
}

/// \brief Check that 1fseB00 has 70 residues
BOOST_FIXTURE_TEST_CASE(prot_1fseB00_has_70_residues, fixture_1a04A02_1fseB00) {
	BOOST_CHECK_EQUAL(prot2.get_length(), 70_z); // 1fseB00
}

/// \brief Check that the 1a04A02/1fseB00 context_sec() scores are as expected
BOOST_FIXTURE_TEST_CASE(context_sec_scores_as_expected_1a04A02_1fseB00, fixture_1a04A02_1fseB00) {
	check_context_sec_scores_as_expected();
}

/// \brief Check that the 1a04A02/1fseB00 residues_have_similar_area_angle_props() scores are as expected
BOOST_FIXTURE_TEST_CASE(residues_have_similar_area_angle_props_1a04A02_1fseB00, fixture_1a04A02_1fseB00) {
	check_residues_have_similar_area_angle_props();
}

BOOST_AUTO_TEST_SUITE_END()

