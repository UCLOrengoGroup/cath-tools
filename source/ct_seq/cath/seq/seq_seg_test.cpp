/// \file
/// \brief The seq_seg test suite

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

#include "cath/common/rapidjson_addenda/to_rapidjson_string.hpp"
#include "cath/seq/seq_seg.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"

using namespace cath::common;
using namespace cath::seq;

using boost::lexical_cast;
using std::invalid_argument;
using std::string;

BOOST_AUTO_TEST_SUITE(seq_seg_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert( get_start_res_index( seq_seg{ 1272, 1363 } ) == 1272, "" );
	static_assert( get_stop_res_index ( seq_seg{ 1272, 1363 } ) == 1363, "" );

	static_assert( get_start_res_index( seq_seg{    0,    1 } ) ==    0, "" );
	static_assert( get_stop_res_index ( seq_seg{    0,    1 } ) ==    1, "" );

	static_assert( get_start_res_index( seq_seg{    0,    0 } ) ==    0, "" );
	static_assert( get_stop_res_index ( seq_seg{    0,    0 } ) ==    0, "" );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(ctor_throws_on_invalid) {
	BOOST_CHECK_THROW( seq_seg( 1, 0 ), invalid_argument );
}

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL( to_string( seq_seg{ 1272, 1320 } ), "seq_seg[1272-1320]" );
}

BOOST_AUTO_TEST_CASE(insetion_operator_works) {
	BOOST_CHECK_EQUAL( lexical_cast<string>( seq_seg{ 1272, 1320 } ), "seq_seg[1272-1320]" );
}

BOOST_AUTO_TEST_CASE(are_overlapping_works) {
	constexpr auto seq_seg_a = seq_seg{ 1266, 1344 };
	constexpr auto seq_seg_b = seq_seg{ 1272, 1320 };
	constexpr auto seq_seg_c = seq_seg{ 1398, 1437 };
	static_assert(   are_overlapping( seq_seg_a, seq_seg_b ), "" );
	static_assert( ! are_overlapping( seq_seg_a, seq_seg_c ), "" );
	static_assert( ! are_overlapping( seq_seg_b, seq_seg_c ), "" );
	static_assert(   are_overlapping( seq_seg_b, seq_seg_a ), "" );
	static_assert( ! are_overlapping( seq_seg_c, seq_seg_a ), "" );
	static_assert( ! are_overlapping( seq_seg_c, seq_seg_b ), "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(length) {
	static_assert( get_length( seq_seg{ 10, 10 } ) ==  1, "" );
	static_assert( get_length( seq_seg{ 10, 20 } ) == 11, "" );
	BOOST_CHECK( true );
}


BOOST_AUTO_TEST_CASE(shorter_and_longer_length) {
	static_assert( shorter_length( seq_seg{ 10, 19 }, seq_seg{ 10, 29 } ) == 10, "" );
	static_assert( longer_length ( seq_seg{ 10, 19 }, seq_seg{ 10, 29 } ) == 20, "" );
	BOOST_CHECK( true );
}


BOOST_AUTO_TEST_CASE(overlap_by_works) {
	static_assert( overlap_by( seq_seg{  2,  3 }, seq_seg{ 4, 7 } ) ==  0, "" );
	static_assert( overlap_by( seq_seg{  3,  4 }, seq_seg{ 4, 7 } ) ==  1, "" );
	static_assert( overlap_by( seq_seg{  4,  5 }, seq_seg{ 4, 7 } ) ==  2, "" );
	static_assert( overlap_by( seq_seg{  5,  6 }, seq_seg{ 4, 7 } ) ==  2, "" );
	static_assert( overlap_by( seq_seg{  6,  7 }, seq_seg{ 4, 7 } ) ==  2, "" );
	static_assert( overlap_by( seq_seg{  7,  8 }, seq_seg{ 4, 7 } ) ==  1, "" );
	static_assert( overlap_by( seq_seg{  8,  9 }, seq_seg{ 4, 7 } ) ==  0, "" );

	static_assert( overlap_by( seq_seg{  4,  7 }, seq_seg{ 4, 7 } ) ==  4, "" );
	static_assert( overlap_by( seq_seg{  2,  9 }, seq_seg{ 4, 7 } ) ==  4, "" );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(fraction_overlap_over_shorter_works) {
	static_assert( fraction_overlap_over_shorter( seq_seg{  2,  3 }, seq_seg{ 4, 7 } ) ==  0.00, "" );
	static_assert( fraction_overlap_over_shorter( seq_seg{  3,  4 }, seq_seg{ 4, 7 } ) ==  0.50, "" );
	static_assert( fraction_overlap_over_shorter( seq_seg{  4,  5 }, seq_seg{ 4, 7 } ) ==  1.00, "" );
	static_assert( fraction_overlap_over_shorter( seq_seg{  5,  6 }, seq_seg{ 4, 7 } ) ==  1.00, "" );
	static_assert( fraction_overlap_over_shorter( seq_seg{  6,  7 }, seq_seg{ 4, 7 } ) ==  1.00, "" );
	static_assert( fraction_overlap_over_shorter( seq_seg{  7,  8 }, seq_seg{ 4, 7 } ) ==  0.50, "" );
	static_assert( fraction_overlap_over_shorter( seq_seg{  8,  9 }, seq_seg{ 4, 7 } ) ==  0.00, "" );

	static_assert( fraction_overlap_over_shorter( seq_seg{  4,  7 }, seq_seg{ 4, 7 } ) ==  1.00, "" );
	static_assert( fraction_overlap_over_shorter( seq_seg{  2,  9 }, seq_seg{ 4, 7 } ) ==  1.00, "" );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(fraction_overlap_over_longer_works) {
	static_assert( fraction_overlap_over_longer( seq_seg{  2,  3 }, seq_seg{ 4, 7 } ) ==  0.00, "" );
	static_assert( fraction_overlap_over_longer( seq_seg{  3,  4 }, seq_seg{ 4, 7 } ) ==  0.25, "" );
	static_assert( fraction_overlap_over_longer( seq_seg{  4,  5 }, seq_seg{ 4, 7 } ) ==  0.50, "" );
	static_assert( fraction_overlap_over_longer( seq_seg{  5,  6 }, seq_seg{ 4, 7 } ) ==  0.50, "" );
	static_assert( fraction_overlap_over_longer( seq_seg{  6,  7 }, seq_seg{ 4, 7 } ) ==  0.50, "" );
	static_assert( fraction_overlap_over_longer( seq_seg{  7,  8 }, seq_seg{ 4, 7 } ) ==  0.25, "" );
	static_assert( fraction_overlap_over_longer( seq_seg{  8,  9 }, seq_seg{ 4, 7 } ) ==  0.00, "" );

	static_assert( fraction_overlap_over_longer( seq_seg{  4,  7 }, seq_seg{ 4, 7 } ) ==  1.00, "" );
	static_assert( fraction_overlap_over_longer( seq_seg{  2,  9 }, seq_seg{ 4, 7 } ) ==  0.50, "" );

	BOOST_CHECK( true );
}


BOOST_AUTO_TEST_SUITE(segments_from_bounds_fn)

BOOST_AUTO_TEST_CASE(builds_correct_for_valid_data) {
	const residx_vec  bounds   = { 0, 0,   1, 1,   20, 25 };
	const seq_seg_vec expected = { seq_seg{ 0, 0 }, seq_seg{ 1, 1 }, seq_seg{ 20, 25 } };
	BOOST_CHECK_EQUAL_RANGES( segments_from_bounds( bounds ), expected );
}

BOOST_AUTO_TEST_CASE(rejects_start_after_stop) {
	const residx_vec  bounds   = { 0, 0,   2, 1,   20, 25 };
	BOOST_CHECK_THROW( segments_from_bounds( bounds ), invalid_argument );
}

BOOST_AUTO_TEST_CASE(rejects_start_not_after_prev_stop_start_after_stop) {
	const residx_vec  bounds   = { 0, 0,   0, 1,   20, 25 };
	BOOST_CHECK_THROW( segments_from_bounds( bounds ), invalid_argument );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_SUITE(write)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_seq_seg) {
	const seq_seg eg_seq_seg = seq_seg{ 10, 20 };
	BOOST_CHECK_EQUAL( to_rapidjson_string< json_style::COMPACT >( eg_seq_seg ), "[10,20]"               );
	BOOST_CHECK_EQUAL( to_rapidjson_string< json_style::PRETTY  >( eg_seq_seg ), "[\n    10,\n    20\n]" );
}

BOOST_AUTO_TEST_CASE(to_json_string_works_for_seq_seg_vec) {
	const seq_seg_vec eg_seq_segs = {
		seq_seg{ 10, 20 },
		seq_seg{ 30, 40 },
	};
	BOOST_CHECK_EQUAL( to_rapidjson_string< json_style::COMPACT >( eg_seq_segs ), "[[10,20],[30,40]]" );
	BOOST_CHECK_EQUAL( to_rapidjson_string< json_style::PRETTY  >( eg_seq_segs ), "[\n    [\n        10,\n        20\n    ],\n    [\n        30,\n        40\n    ]\n]" );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
