/// \file
/// \brief The prc_scores_file test suite

//#include <boost/algorithm/string/predicate.hpp> /// ***** TEMPORARY *****
//#include <boost/filesystem.hpp> /// ***** TEMPORARY *****
//#include <boost/log/trivial.hpp> /// ***** TEMPORARY *****
#include <boost/test/auto_unit_test.hpp>

#include "file/prc_scores_file/prc_scores_entry.hpp"
#include "file/prc_scores_file/prc_scores_file.hpp"

//#include <iomanip> /// ***** TEMPORARY *****
//#include <iostream> /// ***** TEMPORARY *****

using namespace cath::file;
using namespace std;

//using boost::algorithm::ends_with; /// ***** TEMPORARY *****
//using boost::filesystem::directory_iterator; /// ***** TEMPORARY *****
//using boost::filesystem::path; /// ***** TEMPORARY *****

namespace cath {
    namespace test {

        /// \brief The prc_scores_file_test_suite_fixture to assist in testing prc_scores_file
        struct prc_scores_file_test_suite_fixture {
        protected:
            ~prc_scores_file_test_suite_fixture() noexcept = default;

            /// The first 15 results for a PRC results file that has duplicate hits in it
            const string example_file_string = R"(# PRC 1.5.3 (PLAN9, SPACE5, ALL_TRANS), compiled on Mar 28 2007
# Copyright (C) 2002-5 Martin Madera and MRC LMB, Cambridge, UK
# Freely distributed under the GNU General Public License
# 
# Command     : /opt/local/apps/linux/bin/prc/prc -algo viterbi -Emax 1E30 /data1/people/ucbctnl/ticket_914_prcs/prc_layered/1x3zA01.prc /data1/people/ucbctnl/ticket_914_prcs/prc_lib_dir/41ad20f7036355c7802885bafaf4df56.lib /data1/people/ucbctnl/ticket_914_prcs/scanned_results/1x3zA01.results 
# Algorithm   : Viterbi
# Match-match : dot2
# Align mode  : local-local
# Alignments  : none
# Simple stop : 1.5
# Max hits    : 100000
# Max E-value : 1.0e+30
# Started     : Thu Sep 10 10:42:35 2015
# 
# E-value fn  : n_unrel / (1.0 + exp(lambda*reverse + kappa))
#  .. n_unrel : 14338
#  .. lambda  : 2.543026
#  .. kappa   : 0.014341
# 
# hmm1  start1  end1    length1 hit_no  hmm2    start2  end2    length2 simple  reverse  E-value
1x3zA01 2       25      26      1       1x3zA01 2       25      26        15.4    11.4   4.0e-09
1x3zA01 4       20      26      1       1w09A00 71      87      92         7.5     3.8      0.83
1x3zA01 2       25      26      1       1zx2B00 3       26      147        6.0     3.5       2.0
1x3zA01 5       23      26      1       2uxwA04 27      45      130        6.6     3.4       2.7
1x3zA01 2       25      26      1       1gvjA00 88      111     146        5.6     3.2       4.3
1x3zA01 3       12      26      2       1gvjA00 71      80      146        3.8     1.3       457
1x3zA01 17      25      26      2       1zx2B00 119     127     147        3.1     0.7   2.1e+03
1x3zA01 18      22      26      3       1gvjA00 44      48      146        3.1     0.7   2.3e+03
1x3zA01 5       15      26      4       1gvjA00 68      78      146        2.7     0.3   4.8e+03
1x3zA01 5       12      26      3       1zx2B00 132     139     147        2.6     0.1   6.2e+03
1x3zA01 15      20      26      5       1gvjA00 54      59      146        2.5     0.1   6.4e+03
1x3zA01 17      24      26      6       1gvjA00 38      45      146        2.5     0.0   6.9e+03
1x3zA01 18      26      26      2       2uxwA04 32      40      130        3.2    -0.0   7.2e+03
1x3zA01 4       25      26      2       1w09A00 34      55      92         3.6    -0.0   7.2e+03
1x3zA01 2       16      26      3       2uxwA04 46      60      130        3.2    -0.0   7.4e+03
# END)";
            istringstream the_iss{ example_file_string };
        };

    }  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(prc_scores_file_test_suite, cath::test::prc_scores_file_test_suite_fixture)

BOOST_AUTO_TEST_CASE(parses_raw_results_correctly) {
    const auto the_entries = prc_scores_file::parse_prc_scores_file( the_iss );
    BOOST_REQUIRE_EQUAL( the_entries.size(), 15 );
    BOOST_CHECK_EQUAL( the_entries.front().get_name_1(), "1x3zA01" );
	BOOST_CHECK_EQUAL( the_entries.front().get_name_2(), "1x3zA01" );
	BOOST_CHECK_EQUAL( the_entries.back().get_name_1(),  "1x3zA01" );
	BOOST_CHECK_EQUAL( the_entries.back().get_name_2(),  "2uxwA04" );
}

BOOST_AUTO_TEST_CASE(parses_uniqued_results_correctly) {
    const auto the_entries = prc_scores_file::parse_prc_scores_file_fancy( the_iss );
    BOOST_REQUIRE_EQUAL( the_entries.size(), 5 );
    BOOST_CHECK_EQUAL( the_entries.front().get_name_1(), "1x3zA01" );
	BOOST_CHECK_EQUAL( the_entries.front().get_name_2(), "1x3zA01" );
	BOOST_CHECK_EQUAL( the_entries.back().get_name_1(),  "1x3zA01" );
	BOOST_CHECK_EQUAL( the_entries.back().get_name_2(),  "1gvjA00" );

}

// BOOST_AUTO_TEST_CASE(comment_out_this_scan_through_prc_results_directory) {
// 	const auto the_dir = path{ "/cath/mothra-data1/people/ucbctnl/ticket_914_prcs/scanned_results" };
// 	for (const auto& x : directory_iterator( the_dir )) {
// 		if ( ends_with( x.path().filename().string(), ".results.scores" ) ) {
// 			BOOST_LOG_TRIVIAL( warning ) << "About to try to parse " << x.path().string();
// 			const auto the_results = prc_scores_file::parse_prc_scores_file_fancy( x.path() );
// 			cerr << "Parsed " << right << setw( 7 ) << the_results.size() << " unique_results from " << x.path().string() << "\n";
// 		}
// 	}
// }

BOOST_AUTO_TEST_SUITE_END()
