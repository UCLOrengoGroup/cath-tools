/// \file
/// \brief The full_hit_list test suite

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
#include <boost/test/unit_test.hpp>

#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/resolve_hits/full_hit.hpp"
#include "cath/resolve_hits/full_hit_list.hpp"
#include "cath/resolve_hits/full_hit_list_fns.hpp"
#include "cath/resolve_hits/options/spec/crh_score_spec.hpp"
#include "cath/resolve_hits/options/spec/crh_segment_spec.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::rslv;
using namespace ::cath::seq;

using ::boost::lexical_cast;
using ::std::string;

namespace {

		/// \brief The full_hit_list_test_suite_fixture to assist in testing full_hit_list
		struct full_hit_list_test_suite_fixture {
		protected:
			~full_hit_list_test_suite_fixture() noexcept = default;

			/// \brief Make an example full_hit_list for testing
			full_hit_list make_eg_full_hit_list() {
				return full_hit_list{
					{
						full_hit( { seq_seg{ 1266,                        1344 }, }, "label_a", 20.0 ),
						full_hit( { seq_seg{ 1273, 1321 }, seq_seg{ 1399, 1438 }, }, "label_b", 21.0 ),
						full_hit( { seq_seg{ 1101,                        1319 }, }, "label_c", 22.0 ),
						full_hit( { seq_seg{ 1301,                        1321 }, }, "label_d", 23.0 ),
						full_hit( { seq_seg{ 1438,                        1439 }, }, "label_e", 24.0 ),
						full_hit( { seq_seg{ 1272, 1320 }, seq_seg{ 1398, 1437 }, }, "label_f", 25.0 ),
					}
				};
			}

			const full_hit_list eg_full_hit_list = make_eg_full_hit_list();
		};

} // namespace

BOOST_FIXTURE_TEST_SUITE(full_hit_list_test_suite, full_hit_list_test_suite_fixture)

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL( to_string( eg_full_hit_list ), "full_hit_list[6 full_hits]" );
}

BOOST_AUTO_TEST_CASE(insertion_operator_works) {
	BOOST_CHECK_EQUAL( lexical_cast<string>( eg_full_hit_list ), "full_hit_list[6 full_hits]" );
}

BOOST_AUTO_TEST_CASE(get_first_label_works) {
	BOOST_CHECK_EQUAL( front( eg_full_hit_list ).get_label(), "label_a" );
}

BOOST_AUTO_TEST_CASE(get_best_score_works) {
	BOOST_CHECK_EQUAL( *get_best_crh_score( eg_full_hit_list, make_neutral_score_spec() ), 25.0 );
}

BOOST_AUTO_TEST_CASE(get_max_stop_works) {
	BOOST_CHECK_EQUAL( *get_max_stop( eg_full_hit_list ), 1439 );
}

BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_CASE(get_max_stop_works) {
	BOOST_CHECK_EQUAL( to_json_string_with_compact_fullhits( eg_full_hit_list ), R"([
    {"match-id":"label_a","score":20.0,"score-type":"crh-value","boundaries":[[1266,1344]]},
    {"match-id":"label_b","score":21.0,"score-type":"crh-value","boundaries":[[1273,1321],[1399,1438]]},
    {"match-id":"label_c","score":22.0,"score-type":"crh-value","boundaries":[[1101,1319]]},
    {"match-id":"label_d","score":23.0,"score-type":"crh-value","boundaries":[[1301,1321]]},
    {"match-id":"label_e","score":24.0,"score-type":"crh-value","boundaries":[[1438,1439]]},
    {"match-id":"label_f","score":25.0,"score-type":"crh-value","boundaries":[[1272,1320],[1398,1437]]}
])" );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
