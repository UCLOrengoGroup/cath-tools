/// \file
/// \brief The spanning_tree's test suite

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#include "common/boost_addenda/graph/spanning_tree.hpp"
#include "common/pair_insertion_operator.hpp"
#include "common/tuple_insertion_operator.hpp"
#include "common/type_aliases.hpp"

using namespace cath;
using namespace cath::common;

using boost::test_tools::per_element;

BOOST_AUTO_TEST_SUITE(spanning_tree_test_suite)


BOOST_AUTO_TEST_SUITE(make_simple_unweighted_spanning_tree_works)

BOOST_AUTO_TEST_CASE(spanning_tree_for_zero_is_empty) {
	BOOST_TEST ( make_simple_unweighted_spanning_tree( 0 ) == size_size_pair_vec(), per_element{} );
	BOOST_TEST( make_simple_unweighted_spanning_tree( 1 ) == size_size_pair_vec(), per_element{} );
}

BOOST_AUTO_TEST_CASE(spanning_tree_for_one_is_empty) {
	BOOST_TEST( make_simple_unweighted_spanning_tree( 0 ) == size_size_pair_vec(), per_element{} );
	BOOST_TEST( make_simple_unweighted_spanning_tree( 1 ) == size_size_pair_vec(), per_element{} );
}

BOOST_AUTO_TEST_CASE(spanning_tree_for_two_is_zero_one) {
	const size_size_pair_vec expected{ { 0, 1 } };
	BOOST_TEST( make_simple_unweighted_spanning_tree( 2 ) == expected, per_element{} );
}

BOOST_AUTO_TEST_CASE(spanning_tree_for_three_is_zero_one_and_one_two) {
	const size_size_pair_vec expected{ { 0, 1 }, { 1, 2 } };
	BOOST_TEST( make_simple_unweighted_spanning_tree( 3 ) == expected, per_element{} );
}

BOOST_AUTO_TEST_SUITE_END()


/// \brief Check that the correct spanning tree is generated from real data from 3.90.400.10
BOOST_AUTO_TEST_CASE(spanning_tree_for_3_90_400_10) {
	const size_size_doub_tpl_vec edges_and_scores = {
		size_size_doub_tpl{ 0, 1, 85.40 }, // should be 3rd in max spanning tree
		size_size_doub_tpl{ 0, 2, 86.25 },
		size_size_doub_tpl{ 0, 3, 87.96 }, // should be 2nd in max spanning tree
		size_size_doub_tpl{ 1, 2, 85.21 },
		size_size_doub_tpl{ 1, 3, 84.20 },
		size_size_doub_tpl{ 2, 3, 88.34 }, // should be 1st in max spanning tree
	};
	constexpr size_t num_items = 4;

	const size_size_doub_tpl_vec expected_max = {
		size_size_doub_tpl{ 2, 3, 88.34 },
		size_size_doub_tpl{ 0, 3, 87.96 },
		size_size_doub_tpl{ 0, 1, 85.40 },
	};
	BOOST_TEST( calc_max_spanning_tree( edges_and_scores, num_items ) == expected_max, per_element{} );

	const size_size_doub_tpl_vec expected_min = {
		size_size_doub_tpl{ 1, 3, 84.20 },
		size_size_doub_tpl{ 1, 2, 85.21 },
		size_size_doub_tpl{ 0, 1, 85.40 },
	};
	BOOST_TEST( calc_min_spanning_tree( edges_and_scores, num_items ) == expected_min, per_element{} );
}

BOOST_AUTO_TEST_CASE(gets_simple_problem_correct) {
	const size_size_doub_tpl_vec edges_and_scores = {
		size_size_doub_tpl{ 0, 1, -4.0 },
		size_size_doub_tpl{ 0, 2, 80.0 },
		size_size_doub_tpl{ 1, 2, 72.0 },
		size_size_doub_tpl{ 1, 3, 10.0 },
		size_size_doub_tpl{ 2, 3, 61.0 },
	};
	constexpr size_t num_items = 4;

	// Construct the got results and the expected results and then compare
	const size_size_doub_tpl_vec expected_max = {
		size_size_doub_tpl{ 0, 2, 80.0 },
		size_size_doub_tpl{ 1, 2, 72.0 },
		size_size_doub_tpl{ 2, 3, 61.0 }
	};
	BOOST_TEST( calc_max_spanning_tree( edges_and_scores, num_items ) == expected_max, per_element{} );

	const size_size_doub_tpl_vec expected_min = {
		size_size_doub_tpl{ 0, 1, -4.0 },
		size_size_doub_tpl{ 1, 3, 10.0 },
		size_size_doub_tpl{ 2, 3, 61.0 }
	};
	BOOST_TEST( calc_min_spanning_tree( edges_and_scores, num_items ) == expected_min, per_element{} );
}

/// \dot
/// digraph example {
///     graph [ ratio=auto ]
///     node  [ shape=circle, fontname=Helvetica, fontsize=10, style=filled, color=indigo, fontcolor=white ];
///     edge  [ fontname=Helvetica, fontsize=10, arrowhead=open, dir="both" ];
///
///      0 ->  5 [label =  "3" ];
///      1 ->  6 [label =  "7" ];
///      2 ->  6 [label =  "7" ];
///      3 ->  6 [label = "10" ];
///      3 ->  9 [label =  "0" ];
///      4 ->  7 [label =  "1" ];
///      4 ->  8 [label =  "9" ];
///      4 -> 12 [label =  "9" ];
///      5 ->  7 [label =  "8" ];
///      7 ->  9 [label =  "6" ];
///      7 -> 10 [label =  "2" ];
///     10 -> 11 [label =  "5" ];
///     13 ->  9 [label =  "4" ];
/// }
/// \enddot
BOOST_AUTO_TEST_CASE(test_a) {
	const size_size_doub_tpl_vec input = {
		{  0,  5,  3.0 },
		{  1,  6,  7.0 },
		{  2,  6,  7.0 },
		{  3,  6, 10.0 },
		{  3,  9,  0.0 },
		{  4,  7,  1.0 },
		{  4,  8,  9.0 },
		{  4, 12,  9.0 },
		{  5,  7,  8.0 },
		{  7,  9,  6.0 },
		{  7, 10,  2.0 },
		{ 10, 11,  5.0 },
		{ 13,  9,  4.0 },
	};

	const size_size_doub_tpl_vec expected = {
		{  0,  5,  3.0 },
		{  5,  7,  8.0 },
		{  7,  9,  6.0 },
		{ 13,  9,  4.0 },
		{  7, 10,  2.0 },
		{ 10, 11,  5.0 },
		{  4,  7,  1.0 },
		{  4, 12,  9.0 },
		{  4,  8,  9.0 },
		{  3,  9,  0.0 },
		{  3,  6, 10.0 },
		{  2,  6,  7.0 },
		{  1,  6,  7.0 },
	};
	BOOST_TEST( order_spanning_tree_from_start( input, 0 ) == expected, per_element{} );
}

BOOST_AUTO_TEST_SUITE_END()
