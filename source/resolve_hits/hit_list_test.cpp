/// \file
/// \brief The hit_list test suite

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

#include "resolve_hits/hit_list.h"

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::rslv;
using namespace cath::test;

using boost::lexical_cast;
using std::string;

BOOST_TEST_DONT_PRINT_LOG_VALUE( hit_list::const_iterator )

namespace cath {
	namespace test {

		/// \brief The hit_list_test_suite_fixture to assist in testing hit_list
		struct hit_list_test_suite_fixture {
		protected:
			~hit_list_test_suite_fixture() noexcept = default;

			/// \brief Make an example hit_list for testing
			hit_list make_eg_hit_list() {
				/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
				return hit_list{
					{
						make_hit_from_res_indices( { { 1266,                 1344 }, }, 1.0, 0 ),
						make_hit_from_res_indices( { { 1273, 1321 }, { 1399, 1438 }, }, 1.0, 1 ),
						make_hit_from_res_indices( { { 1101,                 1319 }, }, 1.0, 2 ),
						make_hit_from_res_indices( { { 1301,                 1321 }, }, 1.0, 3 ),
						make_hit_from_res_indices( { { 1438,                 1439 }, }, 1.0, 4 ),
						make_hit_from_res_indices( { { 1272, 1320 }, { 1398, 1437 }, }, 1.0, 5 ),
					},
					{
						"label_a",
						"label_b",
						"label_c",
						"label_d",
						"label_e",
						"label_f"
					}
				};
			}

			const hit_list eg_hit_list = make_eg_hit_list();
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(hit_list_test_suite, hit_list_test_suite_fixture)

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL( to_string( eg_hit_list ), "hit_list[6hits]" );
}

BOOST_AUTO_TEST_CASE(insertion_operator_works) {
	BOOST_CHECK_EQUAL( lexical_cast<string>( eg_hit_list ), "hit_list[6hits]" );
}

BOOST_AUTO_TEST_CASE(get_max_stop_works) {
	BOOST_CHECK_EQUAL( *get_max_stop( eg_hit_list ), 1439 );
}

BOOST_AUTO_TEST_CASE(find_first_hit_stopping_at_or_after_works) {
	BOOST_CHECK_EQUAL( get_stop_res_index( *find_first_hit_stopping_at_or_after( eg_hit_list, arrow_after_res( 1320 ) ) ), 1321 );
	BOOST_CHECK_EQUAL( get_stop_res_index( *find_first_hit_stopping_after      ( eg_hit_list, arrow_after_res( 1319 ) ) ), 1321 );

	BOOST_CHECK_EQUAL( get_stop_res_index( *find_first_hit_stopping_at_or_after( eg_hit_list, arrow_after_res( 1330 ) ) ), 1344 );
	BOOST_CHECK_EQUAL( get_stop_res_index( *find_first_hit_stopping_after      ( eg_hit_list, arrow_after_res( 1330 ) ) ), 1344 );

	BOOST_CHECK_EQUAL( get_stop_res_index( *find_first_hit_stopping_at_or_after( eg_hit_list, arrow_after_res( 1438 ) ) ), 1438 );
	BOOST_CHECK_EQUAL( get_stop_res_index( *find_first_hit_stopping_after      ( eg_hit_list, arrow_after_res( 1437 ) ) ), 1438 );

	BOOST_CHECK_EQUAL( get_stop_res_index( *find_first_hit_stopping_at_or_after( eg_hit_list, arrow_after_res( 1439 ) ) ), 1439 );
	BOOST_CHECK_EQUAL( get_stop_res_index( *find_first_hit_stopping_after      ( eg_hit_list, arrow_after_res( 1438 ) ) ), 1439 );

	BOOST_CHECK_EQUAL( find_first_hit_stopping_at_or_after( eg_hit_list, arrow_after_res( 1440 ) ), common::cend( eg_hit_list ) );
	BOOST_CHECK_EQUAL( find_first_hit_stopping_after      ( eg_hit_list, arrow_after_res( 1439 ) ), common::cend( eg_hit_list ) );
}

BOOST_AUTO_TEST_SUITE_END()
