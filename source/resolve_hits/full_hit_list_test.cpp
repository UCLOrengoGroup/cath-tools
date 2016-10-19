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
#include <boost/test/auto_unit_test.hpp>

#include "common/boost_addenda/range/front.h"
#include "resolve_hits/full_hit_list.h"

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace cath::test;

using boost::lexical_cast;
using std::string;

BOOST_TEST_DONT_PRINT_LOG_VALUE( full_hit_list::const_iterator )

namespace cath {
	namespace test {

		/// \brief The full_hit_list_test_suite_fixture to assist in testing full_hit_list
		struct full_hit_list_test_suite_fixture {
		protected:
			~full_hit_list_test_suite_fixture() noexcept = default;

			/// \brief Make an example full_hit_list for testing
			full_hit_list make_eg_full_hit_list() {
				/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
				return full_hit_list{
					{
						full_hit( { hit_seg_of_res_idcs( 1266,                                    1344 ), }, "label_a", 20.0 ),
						full_hit( { hit_seg_of_res_idcs( 1273, 1321 ), hit_seg_of_res_idcs( 1399, 1438 ), }, "label_b", 21.0 ),
						full_hit( { hit_seg_of_res_idcs( 1101,                                    1319 ), }, "label_c", 22.0 ),
						full_hit( { hit_seg_of_res_idcs( 1301,                                    1321 ), }, "label_d", 23.0 ),
						full_hit( { hit_seg_of_res_idcs( 1438,                                    1439 ), }, "label_e", 24.0 ),
						full_hit( { hit_seg_of_res_idcs( 1272, 1320 ), hit_seg_of_res_idcs( 1398, 1437 ), }, "label_f", 25.0 ),
					}
				};
			}

			const full_hit_list eg_full_hit_list = make_eg_full_hit_list();
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(full_hit_list_test_suite, full_hit_list_test_suite_fixture)

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL( to_string( eg_full_hit_list ), "full_hit_list[6full_hits]" );
}

BOOST_AUTO_TEST_CASE(insertion_operator_works) {
	BOOST_CHECK_EQUAL( lexical_cast<string>( eg_full_hit_list ), "full_hit_list[6full_hits]" );
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

BOOST_AUTO_TEST_SUITE_END()
