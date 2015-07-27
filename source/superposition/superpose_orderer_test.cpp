/// \file


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

#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "common/boost_check_no_throw_diag.h"
#include "common/pair_insertion_operator.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "superposition/superpose_orderer.h"

using namespace cath;
using namespace cath::common;
using namespace cath::sup;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The superpose_orderer_test_suite_fixture to assist in testing superpose_orderer
		struct superpose_orderer_test_suite_fixture {
		protected:
			~superpose_orderer_test_suite_fixture() noexcept = default;
		};
	}
}

BOOST_FIXTURE_TEST_SUITE(superpose_orderer_test_suite, cath::test::superpose_orderer_test_suite_fixture)

/// \brief Check that the static private method half_matrix_index_of_indices works for the first few rows and columns
BOOST_AUTO_TEST_CASE(half_matrix_indexing) {
	BOOST_CHECK_EQUAL(  0_z, superpose_orderer::half_matrix_index_of_indices(1, 0) );

	BOOST_CHECK_EQUAL(  1_z, superpose_orderer::half_matrix_index_of_indices(2, 0) );
	BOOST_CHECK_EQUAL(  2_z, superpose_orderer::half_matrix_index_of_indices(2, 1) );

	BOOST_CHECK_EQUAL(  3_z, superpose_orderer::half_matrix_index_of_indices(3, 0) );
	BOOST_CHECK_EQUAL(  4_z, superpose_orderer::half_matrix_index_of_indices(3, 1) );
	BOOST_CHECK_EQUAL(  5_z, superpose_orderer::half_matrix_index_of_indices(3, 2) );

	BOOST_CHECK_EQUAL(  6_z, superpose_orderer::half_matrix_index_of_indices(4, 0) );
	BOOST_CHECK_EQUAL(  7_z, superpose_orderer::half_matrix_index_of_indices(4, 1) );
	BOOST_CHECK_EQUAL(  8_z, superpose_orderer::half_matrix_index_of_indices(4, 2) );
	BOOST_CHECK_EQUAL(  9_z, superpose_orderer::half_matrix_index_of_indices(4, 3) );

	BOOST_CHECK_EQUAL( 10_z, superpose_orderer::half_matrix_index_of_indices(5, 0) );
	BOOST_CHECK_EQUAL( 11_z, superpose_orderer::half_matrix_index_of_indices(5, 1) );
	BOOST_CHECK_EQUAL( 12_z, superpose_orderer::half_matrix_index_of_indices(5, 2) );
	BOOST_CHECK_EQUAL( 13_z, superpose_orderer::half_matrix_index_of_indices(5, 3) );
	BOOST_CHECK_EQUAL( 14_z, superpose_orderer::half_matrix_index_of_indices(5, 4) );
}

/// \brief Check that the static private method half_matrix_index_of_indices works for the first few rows and columns
BOOST_AUTO_TEST_CASE(half_matrix_indexing_throws_for_invalid_indices) {
	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(0, 0), invalid_argument_exception);

	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(0, 1), invalid_argument_exception);
	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(1, 1), invalid_argument_exception);

	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(0, 2), invalid_argument_exception);
	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(1, 2), invalid_argument_exception);
	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(2, 2), invalid_argument_exception);

	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(0, 3), invalid_argument_exception);
	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(1, 3), invalid_argument_exception);
	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(2, 3), invalid_argument_exception);
	BOOST_CHECK_THROW(superpose_orderer::half_matrix_index_of_indices(3, 3), invalid_argument_exception);
}

/// \brief Check the basic functionality of setting and getting scores
BOOST_AUTO_TEST_CASE(set_and_get_score) {
	// Initialise
	superpose_orderer my_superpose_orderer(3);

	// Check that an exception is thrown on attempts to get scores that haven't been set yet
	BOOST_CHECK_THROW(my_superpose_orderer.get_score(1, 0), invalid_argument_exception);
	BOOST_CHECK_THROW(my_superpose_orderer.get_score(2, 0), invalid_argument_exception);
	BOOST_CHECK_THROW(my_superpose_orderer.get_score(2, 1), invalid_argument_exception);

	// Check that scores can be set for valid locations
	BOOST_CHECK_NO_THROW_DIAG(my_superpose_orderer.set_score(1, 0,  0.0));
	BOOST_CHECK_NO_THROW_DIAG(my_superpose_orderer.set_score(2, 0, -1.5));
	BOOST_CHECK_NO_THROW_DIAG(my_superpose_orderer.set_score(2, 1,  1.5));

	// Check that scores can be retrieve from locations where they have been set
	BOOST_CHECK_EQUAL(  0.0, my_superpose_orderer.get_score(1, 0));
	BOOST_CHECK_EQUAL( -1.5, my_superpose_orderer.get_score(2, 0));
	BOOST_CHECK_EQUAL(  1.5, my_superpose_orderer.get_score(2, 1));

	// Check that set_score and get_score both throw an exception for locations out of range
	BOOST_CHECK_THROW(my_superpose_orderer.set_score(3, 0, 0.0), invalid_argument_exception);
	BOOST_CHECK_THROW(my_superpose_orderer.get_score(3, 0     ), invalid_argument_exception);
}

/// \brief Check that get_ordering() throws an exception if it cannot connect everything (size two with no scores)
BOOST_AUTO_TEST_CASE(throws_if_not_connected_two) {
	superpose_orderer my_superpose_orderer_two(2);
	BOOST_CHECK_THROW( get_spanning_tree(my_superpose_orderer_two), invalid_argument_exception );
}

/// \brief Check that get_ordering() throws an exception if it cannot connect everything (size four with two scores)
BOOST_AUTO_TEST_CASE(throws_if_not_connected_four) {
	superpose_orderer my_superpose_orderer_four(4);
	my_superpose_orderer_four.set_score(1, 0, 100.0);
	my_superpose_orderer_four.set_score(3, 2, 100.0);
	BOOST_CHECK_THROW( get_spanning_tree(my_superpose_orderer_four), invalid_argument_exception );
}

/// \brief Check that get_ordering() throws an exception if it cannot connect everything
BOOST_AUTO_TEST_CASE(gets_simple_problem_correct) {
	// Initialise
	superpose_orderer my_superpose_orderer(4);

	// Set some scores
	my_superpose_orderer.set_score(1, 0, -4.0 );
	my_superpose_orderer.set_score(2, 0, 80.0 );
	my_superpose_orderer.set_score(2, 1, 72.0 );
	my_superpose_orderer.set_score(3, 1, 10.0 );
	my_superpose_orderer.set_score(3, 2, 61.0 );

	// Construct the got results and the expected results and then compare
	const size_size_pair_vec got = get_spanning_tree(my_superpose_orderer);
	const size_size_pair_vec expected = {
		{ 0, 2 },
		{ 1, 2 },
		{ 2, 3 }
	};
	BOOST_CHECK_EQUAL_RANGES( expected, got );
}

/// \brief Check that the correct spanning tree is generated from real data from 3.90.400.10
BOOST_AUTO_TEST_CASE(spanning_tree_for_3_90_400_10) {
	const size_size_pair_doub_map data = { { { 0, 1 }, 85.40 },  // should be 3rd in spanning tree
	                                       { { 0, 2 }, 86.25 },
	                                       { { 0, 3 }, 87.96 },  // should be 2nd in spanning tree
	                                       { { 1, 2 }, 85.21 },
	                                       { { 1, 3 }, 84.20 },
	                                       { { 2, 3 }, 88.34 } }; // should be 1st in spanning tree
	size_size_pair_vec got_spanning_tree = get_spanning_tree_ordered_by_desc_score( data );
	const size_size_pair_vec expected = { { 2, 3 },
	                                      { 0, 3 },
	                                      { 0, 1 } };
	BOOST_CHECK_EQUAL_RANGES( expected, got_spanning_tree );
}

BOOST_AUTO_TEST_SUITE_END()

