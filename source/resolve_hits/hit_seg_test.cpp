/// \file
/// \brief The hit_seg test suite

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

#include "resolve_hits/hit_seg.h"

using namespace cath::rslv;

using boost::lexical_cast;
using std::invalid_argument;
using std::string;

BOOST_AUTO_TEST_SUITE(hit_seg_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert( get_start_res_index( hit_seg_of_res_idcs( 1272, 1363 ) ) == 1272, "" );
	static_assert( get_stop_res_index ( hit_seg_of_res_idcs( 1272, 1363 ) ) == 1363, "" );

	static_assert( get_start_res_index( hit_seg_of_res_idcs(    0,    1 ) ) ==    0, "" );
	static_assert( get_stop_res_index ( hit_seg_of_res_idcs(    0,    1 ) ) ==    1, "" );

	static_assert( get_start_res_index( hit_seg_of_res_idcs(    0,    0 ) ) ==    0, "" );
	static_assert( get_stop_res_index ( hit_seg_of_res_idcs(    0,    0 ) ) ==    0, "" );

	BOOST_CHECK_THROW( hit_seg_of_res_idcs( 1, 0 ), invalid_argument );
}

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL( to_string( hit_seg_of_res_idcs( 1272, 1320 ) ), "hit_seg[1272-1320]" );
}

BOOST_AUTO_TEST_CASE(insetion_operator_works) {
	BOOST_CHECK_EQUAL( lexical_cast<string>( hit_seg_of_res_idcs( 1272, 1320 ) ), "hit_seg[1272-1320]" );
}

BOOST_AUTO_TEST_CASE(overlap) {
	constexpr auto hit_seg_a = hit_seg_of_res_idcs( 1266, 1344 );
	constexpr auto hit_seg_b = hit_seg_of_res_idcs( 1272, 1320 );
	constexpr auto hit_seg_c = hit_seg_of_res_idcs( 1398, 1437 );
	static_assert(   hit_segs_overlap( hit_seg_a, hit_seg_b ), "" );
	static_assert( ! hit_segs_overlap( hit_seg_a, hit_seg_c ), "" );
	static_assert( ! hit_segs_overlap( hit_seg_b, hit_seg_c ), "" );
	static_assert(   hit_segs_overlap( hit_seg_b, hit_seg_a ), "" );
	static_assert( ! hit_segs_overlap( hit_seg_c, hit_seg_a ), "" );
	static_assert( ! hit_segs_overlap( hit_seg_c, hit_seg_b ), "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
