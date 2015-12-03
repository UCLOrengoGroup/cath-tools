/// \file
/// \brief The ssap_and_prc test suite

#include <boost/test/auto_unit_test.hpp>

#include "common/boost_check_no_throw_diag.h"
#include "exception/invalid_argument_exception.h"
#include "file/prc_scores_file/prc_scores_entry.h"
#include "file/ssap_scores_file/ssap_scores_entry.h"
#include "score/homcheck_tools/ssap_and_prc.h"

using namespace cath::common;
using namespace cath::file;
using namespace cath::homcheck;

namespace cath {
	namespace test {

		/// \brief The ssap_and_prc_test_suite_fixture to assist in testing ssap_and_prc
		struct ssap_and_prc_test_suite_fixture {
		protected:
			~ssap_and_prc_test_suite_fixture() noexcept = default;

			/// \brief An example PRC  result (with IDs that match those of the SSAP result)
			const prc_scores_entry  the_prc     = prc_scores_entry_from_line( "1i4dA00 4       199     201     1       3cazA00 15      209     219       25.3    16.3   1.6e-11" );

			/// \brief An example SSAP result (with IDs that match those of the PRC  result)
			const ssap_scores_entry the_ssap    = ssap_scores_entry_from_line( "1i4dA00  3cazA00  188  202  78.30  169   83   11   4.86" );

			/// \brief An example SSAP result with a different query ID
			const ssap_scores_entry ssap_diff_q = ssap_scores_entry_from_line( "1i4dB00  3cazA00  185  202  74.36  152   75    7   9.08" );

			/// \brief An example SSAP result with a different match ID
			const ssap_scores_entry ssap_diff_m = ssap_scores_entry_from_line( "1i4dA00  3cazB00  188  207  74.85  176   85    8   7.28" );

			/// \brief An example SSAP result with a different query ID and a different match ID
			const ssap_scores_entry ssap_diff_b = ssap_scores_entry_from_line( "1i4dB00  3cazB00  185  207  78.58  178   85   10   4.04" );

			/// \brief The correct magic function value for the_prc and the_ssap
			static constexpr double MAGIC_FUNCTION_VALUE = 89.095880017344072;
		};
	}
}

constexpr double cath::test::ssap_and_prc_test_suite_fixture::MAGIC_FUNCTION_VALUE;

BOOST_FIXTURE_TEST_SUITE(ssap_and_prc_test_suite, cath::test::ssap_and_prc_test_suite_fixture)

BOOST_AUTO_TEST_CASE(ctor_from_matching_pair_does_not_throw) {
	BOOST_CHECK_NO_THROW_DIAG( ssap_and_prc tmp( the_ssap, the_prc ) );
}

BOOST_AUTO_TEST_CASE(ctor_from_mismatching_pair_throws) {
	BOOST_CHECK_THROW( ssap_and_prc tmp( ssap_diff_q, the_prc ), invalid_argument_exception );
	BOOST_CHECK_THROW( ssap_and_prc tmp( ssap_diff_m, the_prc ), invalid_argument_exception );
	BOOST_CHECK_THROW( ssap_and_prc tmp( ssap_diff_b, the_prc ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(getters_work) {
	const ssap_and_prc the_ssap_and_prc{ the_ssap, the_prc };
	BOOST_CHECK_EQUAL( the_ssap_and_prc.get_ssap(),              the_ssap                   );
	BOOST_CHECK_EQUAL( the_ssap_and_prc.get_prc(),                the_prc                   );
	BOOST_CHECK_EQUAL( the_ssap_and_prc.get_query_id(),           the_prc.get_name_1()      );
	BOOST_CHECK_EQUAL( the_ssap_and_prc.get_match_id(),           the_prc.get_name_2()      );

	BOOST_CHECK_EQUAL( get_ssap_length_1   ( the_ssap_and_prc ), the_ssap.get_length_1   () );
	BOOST_CHECK_EQUAL( get_ssap_length_2   ( the_ssap_and_prc ), the_ssap.get_length_2   () );
	BOOST_CHECK_EQUAL( get_ssap_score      ( the_ssap_and_prc ), the_ssap.get_ssap_score () );
	BOOST_CHECK_EQUAL( get_ssap_num_equivs ( the_ssap_and_prc ), the_ssap.get_num_equivs () );
	BOOST_CHECK_EQUAL( get_ssap_overlap_pc ( the_ssap_and_prc ), the_ssap.get_overlap_pc () );
	BOOST_CHECK_EQUAL( get_ssap_seq_id_pc  ( the_ssap_and_prc ), the_ssap.get_seq_id_pc  () );
	BOOST_CHECK_EQUAL( get_ssap_rmsd       ( the_ssap_and_prc ), the_ssap.get_rmsd       () );

	BOOST_CHECK_EQUAL( get_prc_start_1     ( the_ssap_and_prc ),  the_prc.get_start_1    () );
	BOOST_CHECK_EQUAL( get_prc_end_1       ( the_ssap_and_prc ),  the_prc.get_end_1      () );
	BOOST_CHECK_EQUAL( get_prc_length_1    ( the_ssap_and_prc ),  the_prc.get_length_1   () );
	BOOST_CHECK_EQUAL( get_prc_hit_num     ( the_ssap_and_prc ),  the_prc.get_hit_num    () );
	BOOST_CHECK_EQUAL( get_prc_start_2     ( the_ssap_and_prc ),  the_prc.get_start_2    () );
	BOOST_CHECK_EQUAL( get_prc_end_2       ( the_ssap_and_prc ),  the_prc.get_end_2      () );
	BOOST_CHECK_EQUAL( get_prc_length_2    ( the_ssap_and_prc ),  the_prc.get_length_2   () );
	BOOST_CHECK_EQUAL( get_prc_simple      ( the_ssap_and_prc ),  the_prc.get_simple     () );
	BOOST_CHECK_EQUAL( get_prc_reverse     ( the_ssap_and_prc ),  the_prc.get_reverse    () );
	BOOST_CHECK_EQUAL( get_prc_evalue      ( the_ssap_and_prc ),  the_prc.get_evalue     () );
}

BOOST_AUTO_TEST_CASE(magic_function_works) {
	BOOST_CHECK_EQUAL( magic_function( the_ssap, the_prc ),                          89.095880017344072 );
	BOOST_CHECK_EQUAL( ssap_and_prc( the_ssap, the_prc ).get_magic_function_score(), 89.095880017344072 );
}

BOOST_AUTO_TEST_SUITE_END()
