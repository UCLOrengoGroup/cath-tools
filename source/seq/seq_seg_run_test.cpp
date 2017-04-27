/// \file
/// \brief The seq_seg_run test suite

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

#include "seq/seq_seg_run.hpp"

using namespace cath::seq;

BOOST_AUTO_TEST_SUITE(seq_seg_run_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	const seq_seg_run a{ seq_seg_vec{ {  100,  199 }, {  300,  399 } } };
	const seq_seg_run b{ seq_seg_vec{ {  190,  209 }, {  390,  409 } } };
	const seq_seg_run c{ seq_seg_vec{ {  200,  299 }, {  400,  499 } } };
	const seq_seg_run d{ seq_seg_vec{ { 1100, 1199 }, { 1400, 1499 } } };

	BOOST_CHECK      (   are_overlapping               ( a, a )      );
	BOOST_CHECK      (   are_overlapping               ( a, b )      );
	BOOST_CHECK      ( ! are_overlapping               ( a, c )      );
	BOOST_CHECK      ( ! are_overlapping               ( a, d )      );
	BOOST_CHECK      (   are_overlapping               ( b, a )      );
	BOOST_CHECK      (   are_overlapping               ( b, b )      );
	BOOST_CHECK      (   are_overlapping               ( b, c )      );
	BOOST_CHECK      ( ! are_overlapping               ( b, d )      );
	BOOST_CHECK      ( ! are_overlapping               ( c, a )      );
	BOOST_CHECK      (   are_overlapping               ( c, b )      );
	BOOST_CHECK      (   are_overlapping               ( c, c )      );
	BOOST_CHECK      ( ! are_overlapping               ( c, d )      );
	BOOST_CHECK      ( ! are_overlapping               ( d, a )      );
	BOOST_CHECK      ( ! are_overlapping               ( d, b )      );
	BOOST_CHECK      ( ! are_overlapping               ( d, c )      );
	BOOST_CHECK      (   are_overlapping               ( d, d )      );

	BOOST_CHECK_EQUAL(   overlap_by                    ( a, a ), 200 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( a, b ),  20 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( a, c ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( a, d ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( b, a ),  20 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( b, b ),  40 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( b, c ),  20 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( b, d ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( c, a ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( c, b ),  20 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( c, c ), 200 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( c, d ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( d, a ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( d, b ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( d, c ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( d, d ), 200 );

	BOOST_CHECK_EQUAL(   shorter_length                ( a, a ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( a, b ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( a, c ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( a, d ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( b, a ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( b, b ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( b, c ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( b, d ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( c, a ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( c, b ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( c, c ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( c, d ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( d, a ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( d, b ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( d, c ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( d, d ), 200 );

	BOOST_CHECK_EQUAL(   longer_length                 ( a, a ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( a, b ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( a, c ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( a, d ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( b, a ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( b, b ),  40 );
	BOOST_CHECK_EQUAL(   longer_length                 ( b, c ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( b, d ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( c, a ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( c, b ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( c, c ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( c, d ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( d, a ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( d, b ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( d, c ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( d, d ), 200 );

	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( a, a ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( a, b ), 0.5 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( a, c ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( a, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( b, a ), 0.5 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( b, b ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( b, c ), 0.5 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( b, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( c, a ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( c, b ), 0.5 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( c, c ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( c, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( d, a ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( d, b ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( d, c ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( d, d ), 1.0 );

	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( a, a ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( a, b ), 0.1 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( a, c ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( a, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( b, a ), 0.1 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( b, b ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( b, c ), 0.1 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( b, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( c, a ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( c, b ), 0.1 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( c, c ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( c, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( d, a ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( d, b ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( d, c ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( d, d ), 1.0 );
}

BOOST_AUTO_TEST_CASE(problem_case) {
	const seq_seg_run a{ seq_seg_vec{ {  0, 1 }, { 4, 6 } } };
	const seq_seg_run b{ seq_seg_vec{ {  3, 5 }           } };
	BOOST_CHECK( are_overlapping( a, b ) );
}

BOOST_AUTO_TEST_SUITE_END()
