/// \file
/// \brief The rod_cache test suite

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

#include "common/boost_addenda/range/indices.hpp"
#include "common/boost_addenda/test/boost_check_no_throw_diag.hpp"
#include "common/pair_insertion_operator.hpp"
#include "common/size_t_literal.hpp"
#include "structure/view_cache/detail/plate/rod_cache.hpp"

#include <boost/graph/adjacency_list.hpp>

using namespace cath;
using namespace cath::common;
using namespace cath::index::detail;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The rod_cache_test_suite_fixture to assist in testing rod_cache
		struct rod_cache_test_suite_fixture {
		protected:
			~rod_cache_test_suite_fixture() noexcept = default;
			
		public:

			void check_throw_for_indices(const size_t &,
			                             const size_t &,
			                             const size_t &,
			                             const size_t &) const;
			void check_throw_for_rod_and_notch(const size_t &,
			                                   const size_t &,
			                                   const size_t &,
			                                   const size_t &) const;

		};

	}  // namespace test
}  // namespace cath

/// \brief Check that the code which converts indices to rod and/or notch throws as expected for the specified indices/sizes
///
/// There are six tests:
///  * x2 for do check and don't check
///  * x3 for each of the three functions that get rod and/or notch from indices
void cath::test::rod_cache_test_suite_fixture::check_throw_for_indices(const size_t &arg_index_a, ///< The first index to check
                                                                       const size_t &arg_index_b, ///< The second index to check
                                                                       const size_t &arg_size_a,  ///< The first size of the matrix
                                                                       const size_t &arg_size_b   ///< The second size of the matrix
                                                                       ) const {
	BOOST_CHECK_THROW        ( get_rod_of_indices          ( arg_index_a, arg_index_b, arg_size_a, arg_size_b, true  ), invalid_argument_exception );
	BOOST_CHECK_NO_THROW_DIAG( get_rod_of_indices          ( arg_index_a, arg_index_b, arg_size_a, arg_size_b, false )                             );
	BOOST_CHECK_THROW        ( get_notch_of_indices        ( arg_index_a, arg_index_b, arg_size_a, arg_size_b, true  ), invalid_argument_exception );
	BOOST_CHECK_NO_THROW_DIAG( get_notch_of_indices        ( arg_index_a, arg_index_b, arg_size_a, arg_size_b, false )                             );
	BOOST_CHECK_THROW        ( get_rod_and_notch_of_indices( arg_index_a, arg_index_b, arg_size_a, arg_size_b, true  ), invalid_argument_exception );
	BOOST_CHECK_NO_THROW_DIAG( get_rod_and_notch_of_indices( arg_index_a, arg_index_b, arg_size_a, arg_size_b, false )                             );
}

/// \brief Check that the code which converts rod and notch to index_a and/or index_b throws as expected for the specified indices/sizes
///
/// There are six tests:
///  * x2 for do check and don't check
///  * x3 for each of the three functions that get index_a and/or index_b from rod and notch
void cath::test::rod_cache_test_suite_fixture::check_throw_for_rod_and_notch(const size_t &arg_rod,    ///< The rod value to check
                                                                             const size_t &arg_notch,  ///< The notch value to check
                                                                             const size_t &arg_size_a, ///< The first size of the matrix
                                                                             const size_t &arg_size_b  ///< The second size of the matrix
                                                                             ) const {
	BOOST_CHECK_THROW        ( get_index_a_of_rod_and_notch( arg_rod, arg_notch, arg_size_a, arg_size_b, true  ), invalid_argument_exception );
	BOOST_CHECK_NO_THROW_DIAG( get_index_a_of_rod_and_notch( arg_rod, arg_notch, arg_size_a, arg_size_b, false )                             );
	BOOST_CHECK_THROW        ( get_index_b_of_rod_and_notch( arg_rod, arg_notch, arg_size_a, arg_size_b, true  ), invalid_argument_exception );
	BOOST_CHECK_NO_THROW_DIAG( get_index_b_of_rod_and_notch( arg_rod, arg_notch, arg_size_a, arg_size_b, false )                             );
	BOOST_CHECK_THROW        ( get_indices_of_rod_and_notch( arg_rod, arg_notch, arg_size_a, arg_size_b, true  ), invalid_argument_exception );
	BOOST_CHECK_NO_THROW_DIAG( get_indices_of_rod_and_notch( arg_rod, arg_notch, arg_size_a, arg_size_b, false )                             );
}

/// \brief Test suite to check the rod/cache code
BOOST_FIXTURE_TEST_SUITE(rod_cache_test_suite, cath::test::rod_cache_test_suite_fixture)

/// \brief Check the conversions from indices to rod/cache for a simple 5x4 example
///
/// From indices:
///
///     (0,0)  (0,1)  (0,2)  (0,3)
///     (1,0)  (1,1)  (1,2)  (1,3)
///     (2,0)  (2,1)  (2,2)  (2,3)
///     (3,0)  (3,1)  (3,2)  (3,3)
///     (4,0)  (4,1)  (4,2)  (4,3)
///
/// To rod/cache:
///
///     (4,0)  (5,0)  (6,0)  (7,0)
///     (3,0)  (4,1)  (5,1)  (6,1)
///     (2,0)  (3,1)  (4,2)  (5,2)
///     (1,0)  (2,1)  (3,2)  (4,3)
///     (0,0)  (1,1)  (2,2)  (3,3)
BOOST_AUTO_TEST_CASE(simple_5_by_4_rod_and_notch_of_indices) {
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 0_z, 0_z, 5_z, 4_z ), make_pair( 4_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 0_z, 1_z, 5_z, 4_z ), make_pair( 5_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 0_z, 2_z, 5_z, 4_z ), make_pair( 6_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 0_z, 3_z, 5_z, 4_z ), make_pair( 7_z, 0_z ) );

	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 1_z, 0_z, 5_z, 4_z ), make_pair( 3_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 1_z, 1_z, 5_z, 4_z ), make_pair( 4_z, 1_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 1_z, 2_z, 5_z, 4_z ), make_pair( 5_z, 1_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 1_z, 3_z, 5_z, 4_z ), make_pair( 6_z, 1_z ) );

	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 2_z, 0_z, 5_z, 4_z ), make_pair( 2_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 2_z, 1_z, 5_z, 4_z ), make_pair( 3_z, 1_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 2_z, 2_z, 5_z, 4_z ), make_pair( 4_z, 2_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 2_z, 3_z, 5_z, 4_z ), make_pair( 5_z, 2_z ) );

	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 3_z, 0_z, 5_z, 4_z ), make_pair( 1_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 3_z, 1_z, 5_z, 4_z ), make_pair( 2_z, 1_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 3_z, 2_z, 5_z, 4_z ), make_pair( 3_z, 2_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 3_z, 3_z, 5_z, 4_z ), make_pair( 4_z, 3_z ) );

	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 4_z, 0_z, 5_z, 4_z ), make_pair( 0_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 4_z, 1_z, 5_z, 4_z ), make_pair( 1_z, 1_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 4_z, 2_z, 5_z, 4_z ), make_pair( 2_z, 2_z ) );
	BOOST_CHECK_EQUAL( get_rod_and_notch_of_indices( 4_z, 3_z, 5_z, 4_z ), make_pair( 3_z, 3_z ) );
}

/// \brief Check the conversions from rod/cache to indices and for a simple 5x4 example
///
/// From rod/cache:
///
///     (4,0)  (5,0)  (6,0)  (7,0)
///     (3,0)  (4,1)  (5,1)  (6,1)
///     (2,0)  (3,1)  (4,2)  (5,2)
///     (1,0)  (2,1)  (3,2)  (4,3)
///     (0,0)  (1,1)  (2,2)  (3,3)
///
/// To indices:
///
///     (0,0)  (0,1)  (0,2)  (0,3)
///     (1,0)  (1,1)  (1,2)  (1,3)
///     (2,0)  (2,1)  (2,2)  (2,3)
///     (3,0)  (3,1)  (3,2)  (3,3)
///     (4,0)  (4,1)  (4,2)  (4,3)
BOOST_AUTO_TEST_CASE(simple_5_by_4_indices_of_rod_and_notch) {
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 0_z, 0_z, 5_z, 4_z ), make_pair( 4_z, 0_z ) );

	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 1_z, 0_z, 5_z, 4_z ), make_pair( 3_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 1_z, 1_z, 5_z, 4_z ), make_pair( 4_z, 1_z ) );

	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 2_z, 0_z, 5_z, 4_z ), make_pair( 2_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 2_z, 1_z, 5_z, 4_z ), make_pair( 3_z, 1_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 2_z, 2_z, 5_z, 4_z ), make_pair( 4_z, 2_z ) );

	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 3_z, 0_z, 5_z, 4_z ), make_pair( 1_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 3_z, 1_z, 5_z, 4_z ), make_pair( 2_z, 1_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 3_z, 2_z, 5_z, 4_z ), make_pair( 3_z, 2_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 3_z, 3_z, 5_z, 4_z ), make_pair( 4_z, 3_z ) );

	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 4_z, 0_z, 5_z, 4_z ), make_pair( 0_z, 0_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 4_z, 1_z, 5_z, 4_z ), make_pair( 1_z, 1_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 4_z, 2_z, 5_z, 4_z ), make_pair( 2_z, 2_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 4_z, 3_z, 5_z, 4_z ), make_pair( 3_z, 3_z ) );

	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 5_z, 0_z, 5_z, 4_z ), make_pair( 0_z, 1_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 5_z, 1_z, 5_z, 4_z ), make_pair( 1_z, 2_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 5_z, 2_z, 5_z, 4_z ), make_pair( 2_z, 3_z ) );

	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 6_z, 0_z, 5_z, 4_z ), make_pair( 0_z, 2_z ) );
	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 6_z, 1_z, 5_z, 4_z ), make_pair( 1_z, 3_z ) );

	BOOST_CHECK_EQUAL( get_indices_of_rod_and_notch( 7_z, 0_z, 5_z, 4_z ), make_pair( 0_z, 3_z ) );
}

/// \brief Check some specific out-of-range throws for conversions between indices and rod/notch on a 5x4 matrix
BOOST_AUTO_TEST_CASE(throws_for_5_by_4) {
	check_throw_for_indices( 0_z, 4_z, 5_z, 4_z );
	check_throw_for_indices( 1_z, 4_z, 5_z, 4_z );
	check_throw_for_indices( 2_z, 4_z, 5_z, 4_z );
	check_throw_for_indices( 3_z, 4_z, 5_z, 4_z );
	check_throw_for_indices( 4_z, 4_z, 5_z, 4_z );
	check_throw_for_indices( 5_z, 4_z, 5_z, 4_z );
	check_throw_for_indices( 5_z, 3_z, 5_z, 4_z );
	check_throw_for_indices( 5_z, 2_z, 5_z, 4_z );
	check_throw_for_indices( 5_z, 1_z, 5_z, 4_z );
	check_throw_for_indices( 5_z, 0_z, 5_z, 4_z );

	check_throw_for_rod_and_notch( 0_z, 1_z, 5_z, 4_z );
	check_throw_for_rod_and_notch( 1_z, 2_z, 5_z, 4_z );
	check_throw_for_rod_and_notch( 2_z, 3_z, 5_z, 4_z );
	check_throw_for_rod_and_notch( 3_z, 4_z, 5_z, 4_z );
	check_throw_for_rod_and_notch( 4_z, 4_z, 5_z, 4_z );
	check_throw_for_rod_and_notch( 5_z, 3_z, 5_z, 4_z );
	check_throw_for_rod_and_notch( 6_z, 2_z, 5_z, 4_z );
	check_throw_for_rod_and_notch( 7_z, 1_z, 5_z, 4_z );
	check_throw_for_rod_and_notch( 8_z, 0_z, 5_z, 4_z );
}

/// \brief Test the exceptions for conversion between indices and rod/notch for small matrices
BOOST_AUTO_TEST_CASE(exceptions_for_small_entries) {
	check_throw_for_indices( 0, 0, 0, 1 );
	check_throw_for_indices( 0, 0, 0, 1 );
	check_throw_for_indices( 0, 0, 0, 1 );

	check_throw_for_indices( 0, 0, 1, 0 );
	check_throw_for_indices( 0, 0, 1, 0 );
	check_throw_for_indices( 0, 0, 1, 0 );

	check_throw_for_indices( 0, 0, 0, 0 );
	check_throw_for_indices( 0, 0, 0, 0 );
	check_throw_for_indices( 0, 0, 0, 0 );

	check_throw_for_indices( 1, 0, 1, 1 );
	check_throw_for_indices( 1, 0, 1, 1 );
	check_throw_for_indices( 1, 0, 1, 1 );

	check_throw_for_indices( 0, 1, 1, 1 );
	check_throw_for_indices( 0, 1, 1, 1 );
	check_throw_for_indices( 0, 1, 1, 1 );

	check_throw_for_indices( 1, 1, 1, 1 );
	check_throw_for_indices( 1, 1, 1, 1 );
	check_throw_for_indices( 1, 1, 1, 1 );
}

/// \brief Perform thorough checks on conversions between indices and rod/notch on all matrices up to 20x20
BOOST_AUTO_TEST_CASE(all_sizes_and_indices) {
	const auto max_size = 20_z;
	for (size_t size_a = 1; size_a < max_size; ++size_a) {
		for (size_t size_b = 1; size_b < max_size; ++size_b) {

			// Check that all entries with index_b of 0 have 0 notch and appropriate rod values
			for (const size_t &index_a : indices( size_a ) ) {
				BOOST_CHECK_EQUAL(
					get_rod_and_notch_of_indices( index_a, 0, size_a, size_b ),
					make_pair( size_a - 1 - index_a, 0_z )
				);
			}

			// Check that all entries with index_b of 0 have 0 notch and appropriate rod values
			for (const size_t &index_b : indices( size_b ) ) {
				BOOST_CHECK_EQUAL(
					get_rod_and_notch_of_indices( 0, index_b, size_a, size_b ),
					make_pair( size_a - 1 + index_b, 0_z )
				);
			}

			for (const size_t &index_a : indices( size_a ) ) {
				for (const size_t &index_b : indices( size_b ) ) {
					// Check that converting to rod/notch and back returns the original results
					const size_size_pair  rod_and_notch = get_rod_and_notch_of_indices( index_a, index_b, size_a, size_b );
					const size_t &rod   = rod_and_notch.first;
					const size_t &notch = rod_and_notch.second;
					BOOST_CHECK_EQUAL(
						get_indices_of_rod_and_notch( rod, notch, size_a, size_b ),
						make_pair( index_a, index_b )
					);

					// If both indices are greater than zero then check that taking one of each index
					// reduces the notch by one on the same rod
					if ( index_a > 0 && index_b > 0 ) {
						BOOST_CHECK_EQUAL(
							get_rod_and_notch_of_indices( index_a - 1, index_b - 1, size_a, size_b ),
							make_pair( rod, notch - 1 )
						);
					}
				}
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
