/// \file
/// \brief The hmmer_scores_file test suite

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

#include <boost/test/unit_test.hpp>

#include <spdlog/spdlog.h>

#include "cath/file/hmmer_scores_file/hmmer_scores_entry.hpp"
#include "cath/file/hmmer_scores_file/hmmer_scores_file.hpp"

using namespace ::cath::file;
using namespace ::std;

namespace cath {
    namespace test {

        /// \brief The hmmer_scores_file_test_suite_fixture to assist in testing hmmer_scores_file
        struct hmmer_scores_file_test_suite_fixture {
        protected:
            ~hmmer_scores_file_test_suite_fixture() noexcept = default;

            /// The first 15 results for a PRC results file that has duplicate hits in it
            const string example_file_string = R"(#                                                                                                            --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                                                 accession  query name               accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#                                         ------------------- ----------     -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
cath|4_0_0|3ixfA00/1-137-i5                                   -          cath|4_0_0|102mA00/0-153 -            6.4e-18   63.3   0.8   7.3e-18   63.1   0.8   1.0   1   0   0   1   1   1   1 -
cath|4_0_0|1jliA00/14-125-i5                                  -          cath|4_0_0|102mA00/0-153 -              0.021   12.8   0.5      0.08   10.9   0.5   1.8   1   1   0   1   1   1   0 -
cath|4_0_0|1f1eA00/4-154-i5                                   -          cath|4_0_0|102mA00/0-153 -              0.051   11.4   0.1      0.06   11.2   0.1   1.2   1   0   0   1   1   1   0 -
cath|4_0_0|1iv8A03/305-389-i5                                 -          cath|4_0_0|102mA00/0-153 -                0.1   10.4   0.3      0.17    9.6   0.3   1.4   1   1   0   1   1   1   0 -
cath|4_0_0|2iazA00/0-111-i5                                   -          cath|4_0_0|102mA00/0-153 -               0.68    8.3   0.0       1.5    7.2   0.0   1.5   1   0   0   1   1   1   0 -
cath|4_0_0|1mo9A02/87-131_197-244_270-311-i5                  -          cath|4_0_0|102mA00/0-153 -                  1    6.9   0.0       1.6    6.3   0.0   1.4   1   0   0   1   1   1   0 -
cath|4_0_0|1eucB01/1-20_112-246-i5                            -          cath|4_0_0|102mA00/0-153 -                1.2    6.6   0.0       1.8    6.1   0.0   1.3   1   1   0   1   1   1   0 -
cath|4_0_0|2h5nC00/10-133-i5                                  -          cath|4_0_0|102mA00/0-153 -                1.6    6.4   0.0       2.6    5.8   0.0   1.4   1   0   0   1   1   1   0 -
cath|4_0_0|2i6hA02/84-179-i3                                  -          cath|4_0_0|102mA00/0-153 -                1.9    6.7   0.1       3.9    5.7   0.1   1.6   1   1   0   1   1   1   0 -
cath|4_0_0|1i6kA02/183-293-i5                                 -          cath|4_0_0|102mA00/0-153 -                2.2    6.7   0.2        13    4.2   0.1   1.9   1   1   1   2   2   0   0 -
cath|4_0_0|1e91A00/1-85-i5                                    -          cath|4_0_0|102mA00/0-153 -                2.2    6.3   0.4         7    4.7   0.2   1.8   2   0   0   2   2   1   0 -
cath|4_0_0|2k0mA00/1-104-i5                                   -          cath|4_0_0|102mA00/0-153 -                2.6    6.5   0.5       5.6    5.4   0.1   1.7   1   1   1   2   2   1   0 -
cath|4_0_0|1u7nA00/2-329-i5                                   -          cath|4_0_0|102mA00/0-153 -                2.7    4.8   0.0       5.5    3.8   0.0   1.5   2   0   0   2   2   1   0 -
cath|4_0_0|1vfgB02/135-390-i5                                 -          cath|4_0_0|102mA00/0-153 -                3.2    5.1   0.7       3.2    5.1   0.7   1.1   1   0   0   1   1   1   0 -
cath|4_0_0|1tk1A00/8-265-i5                                   -          cath|4_0_0|102mA00/0-153 -                3.3    4.5   0.0       4.4    4.1   0.0   1.2   1   0   0   1   1   1   0 -
#
# Program:         hmmscan
# Version:         3.1b1 (May 2013)
# Pipeline mode:   SCAN
# Query file:      102mA00.COMBS.fa
# Target file:     incomplete.lib
# Option settings: hmmscan -o /dev/null --tblout 102mA00.a.tbl --noali -E 100000000 --max --cpu 6 incomplete.lib 102mA00.COMBS.fa 
# Current dir:     /cath/homes2/ucbctnl/hmmer3_models_for_assignment_benchmark/pairwise_dataset
# Date:            Tue Aug 11 13:41:06 2015
# [ok]
)";
            istringstream the_iss{ example_file_string };
        };

    }  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(hmmer_scores_file_test_suite, cath::test::hmmer_scores_file_test_suite_fixture)

BOOST_AUTO_TEST_CASE(parses_raw_results_correctly) {
    const auto the_entries = hmmer_scores_file::parse_hmmer_scores_file( the_iss );
    BOOST_REQUIRE_EQUAL( the_entries.size(), 15 );
    BOOST_CHECK_EQUAL( the_entries.front().get_name_1(), "3ixfA00" );
	BOOST_CHECK_EQUAL( the_entries.front().get_name_2(), "102mA00" );
	BOOST_CHECK_EQUAL( the_entries.back().get_name_1(),  "1tk1A00" );
	BOOST_CHECK_EQUAL( the_entries.back().get_name_2(),  "102mA00" );
}

//BOOST_AUTO_TEST_CASE(parses_uniqued_results_correctly) {
//    const auto the_entries = hmmer_scores_file::parse_hmmer_scores_file_fancy( the_iss );
//    BOOST_REQUIRE_EQUAL( the_entries.size(), 5 );
//    BOOST_CHECK_EQUAL( the_entries.front().get_name_1(), "1x3zA01" );
//	BOOST_CHECK_EQUAL( the_entries.front().get_name_2(), "1x3zA01" );
//	BOOST_CHECK_EQUAL( the_entries.back().get_name_1(),  "1x3zA01" );
//	BOOST_CHECK_EQUAL( the_entries.back().get_name_2(),  "1gvjA00" );
//
//}

// BOOST_AUTO_TEST_CASE(comment_out_this_scan_through_hmmer_results_directory) {
// 	const auto the_dir = path{ "/cath/mothra-data1/people/ucbctnl/ticket_914_prcs/scanned_results" };
// 	for (const auto& x : directory_iterator( the_dir )) {
// 		if ( ends_with( x.path().filename().string(), ".results.scores" ) ) {
// 			::spdlog::warn( "About to try to parse {}", x.path().string() );
// 			const auto the_results = hmmer_scores_file::parse_hmmer_scores_file_fancy( x.path() );
// 			cerr << "Parsed " << right << setw( 7 ) << the_results.size() << " unique_results from " << x.path().string() << "\n";
// 		}
// 	}
// }

BOOST_AUTO_TEST_SUITE_END()
