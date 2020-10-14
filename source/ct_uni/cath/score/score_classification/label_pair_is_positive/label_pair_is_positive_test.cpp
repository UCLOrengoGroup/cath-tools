/// \file
/// \brief The label_pair_is_positive test suite

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

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/score/score_classification/label_pair_is_positive/label_pair_is_positive.hpp"
//#include "cath/test/global_test_constants.hpp"

using namespace cath::common;
using namespace cath::score;
using namespace std;

namespace cath {
    namespace test {

        /// \brief The label_pair_is_positive_test_suite_fixture to assist in testing label_pair_is_positive
        struct label_pair_is_positive_test_suite_fixture {
        protected:
            ~label_pair_is_positive_test_suite_fixture() noexcept = default;

            const string EXAMPLE_INPUT = R"(1cukA01 1bvsA01 1
1cukA02 1bvsA02 0
1cukA03 1bvsA03 1
)";
            istringstream example_istream{ EXAMPLE_INPUT};
        };

    }  // namespace test
}  // namespace cath


BOOST_FIXTURE_TEST_SUITE(label_pair_is_positive_test_suite, cath::test::label_pair_is_positive_test_suite_fixture)

BOOST_AUTO_TEST_CASE(parses_correctly) {
	const auto eg_label_pair_is_positive = make_label_pair_is_positive( example_istream );
	BOOST_CHECK_EQUAL( eg_label_pair_is_positive.is_positive( "1cukA01", "1bvsA01" ),  true );
	BOOST_CHECK_EQUAL( eg_label_pair_is_positive.is_positive( "1cukA02", "1bvsA02" ), false );
	BOOST_CHECK_EQUAL( eg_label_pair_is_positive.is_positive( "1cukA03", "1bvsA03" ),  true );
	BOOST_CHECK_THROW( eg_label_pair_is_positive.is_positive( "not", "there" ), invalid_argument_exception );
}

BOOST_AUTO_TEST_SUITE_END()
