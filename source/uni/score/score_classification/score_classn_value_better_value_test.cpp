/// \file
/// \brief The score_classn_value_better_value test suite

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

#include "score/score_classification/score_classn_value.hpp"
#include "score/score_classification/score_classn_value_better_value.hpp"

using namespace cath::score;

namespace cath {
	namespace test {

		/// \brief The score_classn_value_better_value_test_suite_fixture to assist in testing score_classn_value_better_value
		struct score_classn_value_better_value_test_suite_fixture {
		protected:
			~score_classn_value_better_value_test_suite_fixture() noexcept = default;

			const score_classn_value one_false_a = { 1.0, false, "a"};
			const score_classn_value two_false_a = { 2.0, false, "a"};
			const score_classn_value one_true_b  = { 1.0, true,  "b"};
			const score_classn_value two_true_b  = { 2.0, true,  "b"};
		};
	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(score_classn_value_better_value_test_suite, cath::test::score_classn_value_better_value_test_suite_fixture)

/// \brief Check score_classn_value_better_value acts as expected for higher_is_better (ie higher_is_better of true)
BOOST_AUTO_TEST_CASE(higher_is_better) {
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( one_false_a, one_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( one_false_a, two_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( one_false_a, one_true_b  ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( one_false_a, two_true_b  ), false );

	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( two_false_a, one_false_a ), true  );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( two_false_a, two_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( two_false_a, one_true_b  ), true  );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( two_false_a, two_true_b  ), false );

	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( one_true_b,  one_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( one_true_b,  two_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( one_true_b,  one_true_b  ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( one_true_b,  two_true_b  ), false );

	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( two_true_b,  one_false_a ), true  );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( two_true_b,  two_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( two_true_b,  one_true_b  ), true  );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( true  )( two_true_b,  two_true_b  ), false );
}

/// \brief Check score_classn_value_better_value acts as expected for higher_is_worse (ie higher_is_better of false)
BOOST_AUTO_TEST_CASE(higher_is_worse) {
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( one_false_a, one_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( one_false_a, two_false_a ), true  );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( one_false_a, one_true_b  ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( one_false_a, two_true_b  ), true  );

	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( two_false_a, one_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( two_false_a, two_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( two_false_a, one_true_b  ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( two_false_a, two_true_b  ), false );

	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( one_true_b,  one_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( one_true_b,  two_false_a ), true  );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( one_true_b,  one_true_b  ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( one_true_b,  two_true_b  ), true  );

	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( two_true_b,  one_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( two_true_b,  two_false_a ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( two_true_b,  one_true_b  ), false );
	BOOST_CHECK_EQUAL( score_classn_value_better_value( false )( two_true_b,  two_true_b  ), false );
}

BOOST_AUTO_TEST_SUITE_END()

