/// \file
/// \brief The ssap_scores_file test suite

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

#include <boost/algorithm/string/join.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "common/size_t_literal.h"
#include "common/type_aliases.h"
#include "file/ssap_scores_file/ssap_scores_file.h"
#include "file/ssap_scores_file/ssap_scores_entry.h"

using namespace boost::algorithm;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::algorithm::join;

namespace cath {
	namespace test {

		/// \brief The ssap_scores_file_test_suite_fixture to assist in testing ssap_scores_file
		struct ssap_scores_file_test_suite_fixture {
		protected:
			~ssap_scores_file_test_suite_fixture() noexcept = default;

			string DEFAULT_SSAP_SCORE_STRING = R"(# This is a comment line
1n0hA01  1ni4B01  188  186  81.01  162   86   11   3.86
1n0hA01  1ozhA01  188  174  91.92  172   91   33   1.57
1n0hA01  1pvdA01  188  179  87.90  176   93   20   3.12
1n0hA01  1trkA02  188  212  77.94  176   83    6   4.51
1n0hA01  2c3mA01  188  257  79.04  180   70   12   4.05
1n0hA01  2ihtA01  188  179  89.19  172   91   25   2.05
1n0hA01  2o1xA02  188  176  81.57  169   89   13   3.89
1ni4B01  1ozhA01  186  174  80.21  160   86   13   4.00
# This is another comment line
1ni4B01  1pvdA01  186  179  79.17  164   88   11   4.10
1ni4B01  1trkA02  186  212  81.36  167   78   13   2.93
1ni4B01  2c3mA01  186  257  74.93  166   64   11   3.77
1ni4B01  2ihtA01  186  179  79.33  161   86   11   4.20
1ni4B01  2o1xA02  186  176  84.67  169   90   17   2.53
1ozhA01  1pvdA01  174  179  89.22  170   94   20   2.06
1ozhA01  1trkA02  174  212  77.28  162   76    9   4.10
1ozhA01  2c3mA01  174  257  78.83  170   66   12   2.85
1ozhA01  2ihtA01  174  179  90.40  166   92   31   1.53
1ozhA01  2o1xA02  174  176  81.73  161   91   10   3.91
1pvdA01  1trkA02  179  212  77.55  169   79    8   5.23
1pvdA01  2c3mA01  179  257  78.39  177   68    8   3.65
1pvdA01  2ihtA01  179  179  87.46  172   96   21   2.81
1pvdA01  2o1xA02  179  176  80.99  168   93   13   4.15
1trkA02  2c3mA01  212  257  74.85  179   69    8   6.52
1trkA02  2ihtA01  212  179  76.33  166   78   11   4.75
1trkA02  2o1xA02  212  176  84.66  173   81   20   3.04
2c3mA01  2ihtA01  257  179  76.42  176   68    9   4.34
2c3mA01  2o1xA02  257  176  76.82  168   65   11   3.69
2ihtA01  2o1xA02  179  176  80.16  165   92    9   4.10
# Last comment line

)";
			istringstream ssap_scores_iss{ DEFAULT_SSAP_SCORE_STRING };
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(ssap_scores_file_test_suite, cath::test::ssap_scores_file_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(parse_into_ssap_scores_entries) {
	const auto the_entries = ssap_scores_file::parse_ssap_scores_file_simple( ssap_scores_iss );
	BOOST_REQUIRE    ( ! the_entries.empty() );
	BOOST_CHECK_EQUAL( the_entries.size(), 28 );
	BOOST_CHECK_EQUAL( the_entries.front().get_name_1(), "1n0hA01" );
	BOOST_CHECK_EQUAL( the_entries.front().get_name_2(), "1ni4B01" );
	BOOST_CHECK_EQUAL( the_entries.back().get_name_1(),  "2ihtA01" );
	BOOST_CHECK_EQUAL( the_entries.back().get_name_2(),  "2o1xA02" );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	const pair<str_vec, size_size_pair_doub_map> ssap_scores_data = ssap_scores_file::parse_ssap_scores_file(ssap_scores_iss);
	const str_vec                 got_ids    = ssap_scores_data.first;
	const size_size_pair_doub_map got_scores = ssap_scores_data.second;
	const str_vec expected_ids = { "1n0hA01", "1ni4B01", "1ozhA01", "1pvdA01", "1trkA02", "2c3mA01", "2ihtA01", "2o1xA02" };
	BOOST_CHECK_EQUAL_RANGES( expected_ids, got_ids );
	BOOST_CHECK_EQUAL(  28_z, got_scores.size() );
	const size_size_pair zero_and_one = make_pair( 0_z, 1_z );
	BOOST_CHECK_EQUAL( 81.01, got_scores.find( zero_and_one )->second );
}
BOOST_AUTO_TEST_SUITE_END()
