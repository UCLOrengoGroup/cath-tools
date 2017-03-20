/// \file
/// \brief The filter_vs_full_score test suite

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

#include <boost/lexical_cast.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/test_tools.hpp"
#include "score/true_pos_false_neg/classn_outcome.hpp"
#include "structure/view_cache/filter/filter_vs_full_score.hpp"
#include "structure/view_cache/filter/detail/filter_vs_full_score_less.hpp"

using namespace cath::common::test;
using namespace cath::index::filter;
using namespace cath::index::filter::detail;
using namespace cath::score;
using namespace std;

using boost::lexical_cast;

namespace cath {
	namespace test {

		/// \brief The filter_vs_full_score_test_suite_fixture to assist in testing filter_vs_full_score
		struct filter_vs_full_score_test_suite_fixture {
		protected:
			~filter_vs_full_score_test_suite_fixture() noexcept = default;
			
		public:
			/// \brief A standard filter score to use throughout the tests
			static constexpr double STD_FILTER_SCORE = 4.0;

			/// \brief A standard full score to use throughout the tests
			static constexpr double STD_FULL_SCORE   = 6.0;

			/// \brief An example filter_vs_full_score to be used in the tests
			filter_vs_full_score std_score_pair = { STD_FILTER_SCORE, STD_FULL_SCORE };
		};

	}  // namespace test
}  // namespace cath

constexpr double cath::test::filter_vs_full_score_test_suite_fixture::STD_FILTER_SCORE;
constexpr double cath::test::filter_vs_full_score_test_suite_fixture::STD_FULL_SCORE;

/// \brief A test suite to unit test filter_vs_full_score
BOOST_FIXTURE_TEST_SUITE(filter_vs_full_score_test_suite, cath::test::filter_vs_full_score_test_suite_fixture)

/// \brief Check the equality/inequality operators work as expected for filter_vs_full_score
BOOST_AUTO_TEST_CASE(equality) {
	const auto score_pairs = {
		std_score_pair,
		filter_vs_full_score( STD_FILTER_SCORE,       STD_FULL_SCORE + 0.5 ),
		filter_vs_full_score( STD_FILTER_SCORE - 0.5, STD_FULL_SCORE       )
	};
	check_equality_operators_on_diff_vals_range( score_pairs );
	check_equality_operators_on_equal_vals( std_score_pair, filter_vs_full_score( STD_FILTER_SCORE,       STD_FULL_SCORE       ) );
}

/// \brief Check the filter_score_less and full_score_less functors for performing less-than operations
BOOST_AUTO_TEST_CASE(less_than) {
	BOOST_CHECK_EQUAL( filter_score_less()( std_score_pair, std_score_pair                                                       ), false );
	BOOST_CHECK_EQUAL( filter_score_less()( std_score_pair, filter_vs_full_score( STD_FILTER_SCORE - 0.5, STD_FULL_SCORE + 0.5 ) ), false );
	BOOST_CHECK_EQUAL( filter_score_less()( std_score_pair, filter_vs_full_score( STD_FILTER_SCORE + 0.5, STD_FULL_SCORE - 0.5 ) ), true  );

	BOOST_CHECK_EQUAL( full_score_less  ()( std_score_pair, std_score_pair                                                       ), false );
	BOOST_CHECK_EQUAL( full_score_less  ()( std_score_pair, filter_vs_full_score( STD_FILTER_SCORE - 0.5, STD_FULL_SCORE + 0.5 ) ), true  );
	BOOST_CHECK_EQUAL( full_score_less  ()( std_score_pair, filter_vs_full_score( STD_FILTER_SCORE + 0.5, STD_FULL_SCORE - 0.5 ) ), false );
}

/// \brief Check the insertion and extraction operators
BOOST_AUTO_TEST_CASE(insertion_and_extraction) {
	BOOST_CHECK_EQUAL(
		lexical_cast<string>( std_score_pair ),
		"filter_vs_full_score[      4,      6]"
	);
	BOOST_CHECK_EQUAL(
		lexical_cast<filter_vs_full_score>( string( "filter_vs_full_score[      4,      6]" ) ),
		std_score_pair
	);
	BOOST_CHECK_EQUAL(
		lexical_cast<filter_vs_full_score>( lexical_cast<string>( std_score_pair ) ),
		std_score_pair
	);
}

/// \brief Check that assess_real_scores_on_filter_attempt calculates correct classn_outcome values for the nine basic categories of case
///        (filter is lt, eq, gt; full is lt, eq, gt)
BOOST_AUTO_TEST_CASE(assess_filter) {
	BOOST_CHECK_EQUAL( assess_real_scores_on_filter_attempt( filter_vs_full_score( STD_FILTER_SCORE + 0.5, STD_FULL_SCORE + 0.5 ), std_score_pair ), classn_outcome::TRUE_POSITIVE  );
	BOOST_CHECK_EQUAL( assess_real_scores_on_filter_attempt( filter_vs_full_score( STD_FILTER_SCORE + 0.5, STD_FULL_SCORE       ), std_score_pair ), classn_outcome::TRUE_POSITIVE  );
	BOOST_CHECK_EQUAL( assess_real_scores_on_filter_attempt( filter_vs_full_score( STD_FILTER_SCORE + 0.5, STD_FULL_SCORE - 0.5 ), std_score_pair ), classn_outcome::FALSE_POSITIVE );
	BOOST_CHECK_EQUAL( assess_real_scores_on_filter_attempt( filter_vs_full_score( STD_FILTER_SCORE,       STD_FULL_SCORE + 0.5 ), std_score_pair ), classn_outcome::TRUE_POSITIVE  );
	BOOST_CHECK_EQUAL( assess_real_scores_on_filter_attempt( filter_vs_full_score( STD_FILTER_SCORE,       STD_FULL_SCORE       ), std_score_pair ), classn_outcome::TRUE_POSITIVE  );
	BOOST_CHECK_EQUAL( assess_real_scores_on_filter_attempt( filter_vs_full_score( STD_FILTER_SCORE,       STD_FULL_SCORE - 0.5 ), std_score_pair ), classn_outcome::FALSE_POSITIVE );
	BOOST_CHECK_EQUAL( assess_real_scores_on_filter_attempt( filter_vs_full_score( STD_FILTER_SCORE - 0.5, STD_FULL_SCORE + 0.5 ), std_score_pair ), classn_outcome::FALSE_NEGATIVE );
	BOOST_CHECK_EQUAL( assess_real_scores_on_filter_attempt( filter_vs_full_score( STD_FILTER_SCORE - 0.5, STD_FULL_SCORE       ), std_score_pair ), classn_outcome::FALSE_NEGATIVE );
	BOOST_CHECK_EQUAL( assess_real_scores_on_filter_attempt( filter_vs_full_score( STD_FILTER_SCORE - 0.5, STD_FULL_SCORE - 0.5 ), std_score_pair ), classn_outcome::TRUE_NEGATIVE  );
}

BOOST_AUTO_TEST_SUITE_END()
