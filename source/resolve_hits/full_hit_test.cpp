/// \file
/// \brief The full_hit test suite

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

#include "common/rapidjson_addenda/to_rapidjson_string.hpp"
#include "resolve_hits/full_hit.hpp"

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace cath::test;

namespace cath {
	namespace test {

		/// \brief The full_hit_test_suite_fixture to assist in testing full_hit
		struct full_hit_test_suite_fixture {
		protected:
			~full_hit_test_suite_fixture() noexcept = default;

			const full_hit eg_full_hit_a{ { hit_seg_of_res_idcs( 1272, 1363 ) }, "lemur", 1.0 };

			const full_hit eg_full_hit_b{ { hit_seg_of_res_idcs( 1272, 1320 ), hit_seg_of_res_idcs( 1398, 1437 ) }, "pangolin", 1.0 };
		};

	}  // namespace test
}  // namespace cath


BOOST_FIXTURE_TEST_SUITE(full_hit_test_suite, full_hit_test_suite_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( get_start_res_index_of_segment( eg_full_hit_a, 0 ), 1272 );
	BOOST_CHECK_EQUAL( get_stop_res_index_of_segment ( eg_full_hit_a, 0 ), 1363 );
	BOOST_CHECK_EQUAL( to_string( eg_full_hit_a ), R"(full_hit[1272-1363; score: 1; label: "lemur"])" );
}

BOOST_AUTO_TEST_CASE(basic_2) {
	BOOST_CHECK_EQUAL( get_start_res_index_of_segment( eg_full_hit_b, 0 ), 1272 );
	BOOST_CHECK_EQUAL( get_stop_res_index_of_segment ( eg_full_hit_b, 0 ), 1320 );
	BOOST_CHECK_EQUAL( get_start_res_index_of_segment( eg_full_hit_b, 1 ), 1398 );
	BOOST_CHECK_EQUAL( get_stop_res_index_of_segment ( eg_full_hit_b, 1 ), 1437 );
	BOOST_CHECK_EQUAL( to_string( eg_full_hit_b ), R"(full_hit[1272-1320,1398-1437; score: 1; label: "pangolin"])" );
}

BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_CASE(get_max_stop_works) {
	BOOST_CHECK_EQUAL(
		to_rapidjson_string<json_style::COMPACT>( eg_full_hit_a ),
		R"({"match-id":"lemur","score":1.0,"score-type":"crh-value","boundaries":[[1272,1363]]})"
	);
	BOOST_CHECK_EQUAL(
		to_rapidjson_string<json_style::COMPACT>( eg_full_hit_b ),
		R"({"match-id":"pangolin","score":1.0,"score-type":"crh-value","boundaries":[[1272,1320],[1398,1437]]})"
	);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()
