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

#include "resolve_hits/full_hit.hpp"

using namespace cath;
using namespace cath::rslv;

BOOST_AUTO_TEST_SUITE(full_hit_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	const auto the_full_hit = full_hit( { hit_seg_of_res_idcs( 1272, 1363 ) }, "lemur", 1.0 );
	BOOST_CHECK_EQUAL( get_start_res_index_of_segment( the_full_hit, 0 ), 1272 );
	BOOST_CHECK_EQUAL( get_stop_res_index_of_segment ( the_full_hit, 0 ), 1363 );
	BOOST_CHECK_EQUAL( to_string( the_full_hit ), R"(full_hit[1272-1363; score: 1; label: "lemur"])" );
}

BOOST_AUTO_TEST_CASE(basic_2) {
	const auto the_full_hit = full_hit( { hit_seg_of_res_idcs( 1272, 1320 ), hit_seg_of_res_idcs( 1398, 1437 ) }, "pangolin", 1.0 );
	BOOST_CHECK_EQUAL( get_start_res_index_of_segment( the_full_hit, 0 ), 1272 );
	BOOST_CHECK_EQUAL( get_stop_res_index_of_segment ( the_full_hit, 0 ), 1320 );
	BOOST_CHECK_EQUAL( get_start_res_index_of_segment( the_full_hit, 1 ), 1398 );
	BOOST_CHECK_EQUAL( get_stop_res_index_of_segment ( the_full_hit, 1 ), 1437 );
	BOOST_CHECK_EQUAL( to_string( the_full_hit ), R"(full_hit[1272-1320,1398-1437; score: 1; label: "pangolin"])" );
}

BOOST_AUTO_TEST_SUITE_END()
