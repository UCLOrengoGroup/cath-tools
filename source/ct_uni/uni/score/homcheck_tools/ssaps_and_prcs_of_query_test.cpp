/// \file
/// \brief The ssaps_and_prcs_of_query test suite

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

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/one_of.hpp>
#include <boost/test/unit_test.hpp>

#include "common/boost_addenda/log/stringstream_log_sink.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/type_aliases.hpp"
#include "file/prc_scores_file/prc_scores_entry.hpp"
#include "file/prc_scores_file/prc_scores_file.hpp"
#include "file/ssap_scores_file/ssap_scores_entry.hpp"
#include "file/ssap_scores_file/ssap_scores_file.hpp"
#include "score/homcheck_tools/first_result_if.hpp"
#include "score/homcheck_tools/ssap_and_prc.hpp"
#include "score/homcheck_tools/ssaps_and_prcs_of_query.hpp"

using namespace cath::file;
using namespace cath::homcheck;
using namespace std;

using boost::algorithm::all_of;
using boost::algorithm::contains;
using boost::algorithm::one_of;
using cath::common::invalid_argument_exception;
using cath::stringstream_log_sink;

namespace cath {
	namespace test {

		/// \brief The ssaps_and_prcs_of_query_test_suite_fixture to assist in testing ssaps_and_prcs_of_query
		struct ssaps_and_prcs_of_query_test_suite_fixture {
		protected:
			~ssaps_and_prcs_of_query_test_suite_fixture() noexcept = default;

			// SOME NOTES:
			//
			// "/export/people/ucbctnl/ticket_914_data/data_data.txt"
			//
			// psql cathdb_current -A -t -F ' ' -c "SELECT domain_id, SUBPATH( cath_id, 0, 4 ) FROM domain WHERE flow_stage_type = 'ASSIGNED' ORDER BY domain_id;" > superfamily_of_domain.txt
			//
			// Magic function >= 79.7779328254 gives an error rate of 5% for a coverage of 46.959%
			// SVM            >= 3.47554072714 gives an error rate of 5% for a coverage of 47.571%

			/// \brief The query ID used in the other valid test data
			const string query_id = "1fkmA02";

			/// \brief The match IDs use in the other valid test data
			const str_vec matching_match_ids = {
				"1pauB00",
				"1qtnB00",
				"1pyoB00",
				"1iceB00",
				"2fdrA02"
			};

			/// \brief Some raw text SSAP scores data
			const string raw_ssap_data = R"(1fkmA02  1pauB00  121   91  67.58   52   42    5   6.33
1fkmA02  1qtnB00  121   90  67.14   56   46    4   7.95
1fkmA02  1pyoB00  121   98  66.72   55   45    4   9.33
1fkmA02  1iceB00  121   88  65.35   49   40    3   9.31
1fkmA02  2fdrA02  121   64  63.67   62   51    6   8.78)";

			/// \brief Some raw text PRC scores data
			const string raw_prc_data = R"(1fkmA02 18      43      141     1       1mi3A00 240     265     322        8.7     4.8      0.82
1fkmA02 15      43      141     1       2bgsA00 210     238     308        8.4     4.6       1.3
1fkmA02 41      50      141     1       1pyoB00 38      47      105        5.0     0.6   3.0e+03
1fkmA02 75      79      141     1       1l6zA02 71      75      109        4.4     0.6   3.0e+03
1fkmA02 79      102     141     1       1pauB00 43      66      102        5.0     0.6   3.0e+03
1fkmA02 23      35      141     1       2fdrA02 8       20      67         3.6     0.4   4.2e+03
1fkmA02 44      63      141     1       1b79B00 28      47      103        4.7     0.4   4.2e+03
1fkmA02 41      50      141     1       1qtnB00 36      45      95         5.0     0.0   6.1e+03
1fkmA02 60      66      141     1       1h2bB01 57      63      219        3.1     0.0   6.1e+03
1fkmA02 41      51      141     1       1iceB00 32      42      88         4.5    -0.3   7.8e+03)";

			/// \brief A ssap_scores_entry_vec parsed from the above raw SSAP scores text
			const ssap_scores_entry_vec ssap_scores = ssap_scores_file::parse_ssap_scores_file_simple( raw_ssap_data );

			/// \brief A prc_scores_entry_vec parsed from the above raw PRC scores text
			const prc_scores_entry_vec  prc_scores  = prc_scores_file::parse_prc_scores_file_fancy   ( raw_prc_data  );

			/// \brief Grab any logging
			const stringstream_log_sink log_sink{};
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(ssaps_and_prcs_of_query_test_suite, cath::test::ssaps_and_prcs_of_query_test_suite_fixture)

BOOST_AUTO_TEST_CASE(building_from_ssaps_and_prcs_works) {
	const auto the_ssaps_and_prcs = make_ssaps_and_prcs_of_query(
		ssap_scores,
		prc_scores
	);

	BOOST_CHECK_EQUAL( the_ssaps_and_prcs.size(), 5 );
	BOOST_CHECK( all_of( the_ssaps_and_prcs, [&] (const ssap_and_prc &x) { return ( x.get_query_id() == query_id ); } ) );
	for (const string &match_id : matching_match_ids) {
		BOOST_CHECK( one_of( the_ssaps_and_prcs, [&] (const ssap_and_prc &x) { return ( x.get_match_id() == match_id ); } ) );
	}
	BOOST_CHECK( contains( log_sink.str(), "5 unmatched PRC results" ) );
}

BOOST_AUTO_TEST_CASE(ctor_throws_on_conflicting_query_ids) {
	const ssap_and_prc_vec ssap_and_prc_entries = {
		ssap_and_prc{
			ssap_scores_entry_from_line( "1fkmA02  1pauB00  121   91  67.58   52   42    5   6.33" ),
			prc_scores_entry_from_line ( "1fkmA02 79      102     141     1       1pauB00 43      66      102        5.0     0.6   3.0e+03" )
		},
		ssap_and_prc{
			ssap_scores_entry_from_line( "1fkmA03  1qtnB00  121   90  67.14   56   46    4   7.95" ),
			prc_scores_entry_from_line ( "1fkmA03 41      50      141     1       1qtnB00 36      45      95         5.0     0.0   6.1e+03" )
		}
	};
	BOOST_CHECK_THROW( ssaps_and_prcs_of_query tmp{ ssap_and_prc_entries }, invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(factory_throws_on_conflicting_query_ids) {
	const auto ssaps = ssap_scores_file::parse_ssap_scores_file_simple( string( "1fkmA02  1pauB00  121   91  67.58   52   42    5   6.33\n1fkmA03  1qtnB00  121   90  67.14   56   46    4   7.95\n" ) );
	const auto prcs  =  prc_scores_file::parse_prc_scores_file_fancy  ( string( "1fkmA02 79      102     141     1       1pauB00 43      66      102        5.0     0.6   3.0e+03\n1fkmA03 41      50      141     1       1qtnB00 36      45      95         5.0     0.0   6.1e+03\n" ) );
	BOOST_CHECK_THROW( make_ssaps_and_prcs_of_query( ssaps, prcs ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(best_magic_function_works) {
	const auto the_ssaps_and_prcs = make_ssaps_and_prcs_of_query(
		ssap_scores,
		prc_scores
	);
	const auto result = first_result_if(
		the_ssaps_and_prcs,
		[] (const ssap_and_prc &x, const ssap_and_prc &y) {
			// Reverse inequality to put the highest magic_function values to the start
			return x.get_magic_function_score() > y.get_magic_function_score();
		},
		[] (const ssap_and_prc &) { return true; }
	);
	BOOST_REQUIRE( result );
	BOOST_CHECK_EQUAL( result->get().get_match_id(), "1pauB00" );
}

BOOST_AUTO_TEST_CASE(best_magic_function_if_works) {
	const auto the_ssaps_and_prcs = make_ssaps_and_prcs_of_query(
		ssap_scores,
		prc_scores
	);
	const auto result = first_result_if(
		the_ssaps_and_prcs,
		[] (const ssap_and_prc &x, const ssap_and_prc &y) {
			// Reverse inequality to put the highest magic_function values to the start
			return x.get_magic_function_score() > y.get_magic_function_score();
		},
		[] (const ssap_and_prc &x) { return x.get_match_id() == "1iceB00"; }
	);
	BOOST_REQUIRE( result );
	BOOST_CHECK_EQUAL( result->get().get_match_id(), "1iceB00" );
}

BOOST_AUTO_TEST_SUITE_END()
