/// \file
/// \brief The windowed_matrix test suite

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

#include "windowed_matrix.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::std;

using ::boost::irange;
using ::boost::numeric_cast;
using ::boost::range::find;
using ::std::crbegin;
using ::std::crend;

namespace {

	/// \brief A fixture for the windowed_matrix_test_suite
	struct windowed_matrix_test_suite_fixture {
	protected:
		~windowed_matrix_test_suite_fixture() noexcept = default;

	public:

		/// \brief Return the second dimensions size of a test matrix
		///
		/// \pre Matrix cannot be empty (ie first dimension size must be > 0) else an invalid_argument_exception will be thrown
		///
		/// \pre Matrix must have a consistent size of the second dimension else an invalid_argument_exception will be thrown
		///
		/// \returns The size of the second dimension
		static size_t get_second_dimension_size(const bool_deq_vec &prm_test_matrix ///< The test matrix to query
		                                              ) {
			if (prm_test_matrix.empty()) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get second dimension size for empty test matrix"));
			}
			const bool_deq::size_type second_dimension_size = prm_test_matrix.front().size();
			for (const bool_deq &entry : prm_test_matrix) {
				if (second_dimension_size != entry.size()) {
					BOOST_THROW_EXCEPTION(invalid_argument_exception("Test matrix does not have a consistent second dimension size"));
				}
			}
			return second_dimension_size;
		}

		/// \brief Transpose one of the test matrices so that the transpose can be tested too
		///
		/// \returns A transposed_copy of the test matrix
		static bool_deq_vec transpose_test_matrix(const bool_deq_vec &prm_test_matrix ///< The test matrix to transpose
		                                                  ) {
			const size_t second_dimension_size = get_second_dimension_size(prm_test_matrix);
			bool_deq_vec transposed_matrix(second_dimension_size);
			for (const bool_deq &entry : prm_test_matrix) {
				for (const size_t &second_dim_ctr : indices( second_dimension_size ) ) {
					transposed_matrix[second_dim_ctr].push_back(entry[second_dim_ctr]);
				}
			}
			return transposed_matrix;
		}

		/// Return the indices (offset 0) of the first and last trues in the list of bools
		///
		/// \pre The deque<bool> must contain at least one true value
		///
		/// \returns A pair<size_t, size_t> containing the indices of the first and last true values respectively
		static size_size_pair get_indices_of_first_and_last_trues(const bool_deq &prm_entry ///< The list of bools to query
		                                                          ) {
			if ( ! contains( prm_entry, true) ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot find indices of first and last true values because there are no true values"));
			}
			const size_t first_true_idx = numeric_cast<size_t>( distance( cbegin ( prm_entry ), find( prm_entry,            true ) ) );
			const size_t last_true_ridx = numeric_cast<size_t>( distance( crbegin( prm_entry ), find( crbegin( prm_entry ), crend( prm_entry), true ) ) );
			const size_t last_true_idx  = prm_entry.size() - last_true_ridx - 1;
	//		cerr << "First and last indices for ";
	//		for (const bool &value : prm_entry) {
	//			cerr << " " << boolalpha << value;
	//		}
	//		cerr << " are " << first_true_idx << " and " << last_true_idx << endl;
			return make_pair(first_true_idx, last_true_idx);
		}

		/// \brief Check that get_window_start_a_for_b() and get_window_stop_a_for_b() recreate the values in the matrix for the specified window size
		void check_windowed_matrix(const bool_deq_vec &prm_test_matrix, ///< The test matrix to check
		                           const size_t       &prm_window_size  ///< The window size to use
		                           ) const {
			const size_t first_dimension_size  = prm_test_matrix.size();
			const size_t second_dimension_size = get_second_dimension_size(prm_test_matrix);

			for (const size_t &entry_ctr : indices( first_dimension_size ) ) {
				const bool_deq       &entry            = prm_test_matrix[entry_ctr];
				const size_size_pair  first_and_last   = get_indices_of_first_and_last_trues(entry);
				const size_t         &first_true_index = first_and_last.first;
				const size_t         &last_true_index  = first_and_last.second;

				BOOST_CHECK_EQUAL(
					first_true_index + 1,
					get_window_start_a_for_b__offset_1(
						second_dimension_size,
						first_dimension_size,
						prm_window_size,
						entry_ctr + 1
					)
				);
				BOOST_CHECK_EQUAL(
					last_true_index + 1,
					get_window_stop_a_for_b__offset_1(
						second_dimension_size,
						first_dimension_size,
						prm_window_size,
						entry_ctr + 1
					)
				);
				if (first_true_index > 0 && last_true_index + 1 < second_dimension_size) {
					if (first_dimension_size != second_dimension_size) {
						BOOST_CHECK_EQUAL(
							prm_window_size,
							last_true_index - first_true_index + 1
						);
					}
				}
			}
		}

		const bool_deq_vec NORMAL_MATRIX = {
			{ true,  true,  true,  true,  true,  false, false, false, false },
			{ true,  true,  true,  true,  true,  true,  false, false, false },
			{ true,  true,  true,  true,  true,  true,  true,  false, false },
			{ false, true,  true,  true,  true,  true,  true,  true,  false },
			{ false, false, true,  true,  true,  true,  true,  true,  true  },
			{ false, false, false, true,  true,  true,  true,  true,  true  },
			{ false, false, false, false, true,  true,  true,  true,  true  }
		};

		const bool_deq_vec SQUARE_MATRIX = {
			{ true,  true,  false, false, false },
			{ true,  true,  true,  false, false },
			{ false, true,  true,  true,  false },
			{ false, false, true,  true,  true  },
			{ false, false, false, true,  true  }
		};

		const bool_deq_vec NARROWEST_SYM_WINDOW = {
			{ true,  true,  true,  false, false },
			{ false, true,  true,  true,  false },
			{ false, false, true,  true,  true  }
		};

		const bool_deq_vec NARROWEST_ASYM_WINDOW = {
			{ true,  true,  false, false, false },
			{ true,  true,  true,  false, false },
			{ false, true,  true,  true,  false },
			{ false, false, true,  true,  true  }
		};

		const bool_deq_vec FULL = {
			{ true , true , true , true },
			{ true , true , true , true },
			{ true , true , true , true }
		};
	};

} // namespace

BOOST_FIXTURE_TEST_SUITE(windowed_matrix_test_suite, windowed_matrix_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(normal) {
	check_windowed_matrix(                       NORMAL_MATRIX,            7 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(normal_transpose) {
	check_windowed_matrix( transpose_test_matrix(NORMAL_MATRIX),           7 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(square) {
	check_windowed_matrix(                       SQUARE_MATRIX,            3 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(square_transpose) {
	check_windowed_matrix( transpose_test_matrix(SQUARE_MATRIX),           3 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(square_odd) {
	check_windowed_matrix(                       SQUARE_MATRIX,            2 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(square_odd_transpose) {
	check_windowed_matrix( transpose_test_matrix(SQUARE_MATRIX),           2 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(narrowest_sym) {
	check_windowed_matrix(                       NARROWEST_SYM_WINDOW,     3 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(narrowest_sym_transpose) {
	check_windowed_matrix( transpose_test_matrix(NARROWEST_SYM_WINDOW),    3 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(narrowest_asym) {
	check_windowed_matrix(                       NARROWEST_ASYM_WINDOW,    3 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(narrowest_asym_transpose) {
	check_windowed_matrix( transpose_test_matrix(NARROWEST_ASYM_WINDOW),   3 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(full) {
	check_windowed_matrix(                       FULL,                   100 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(full_transpose) {
	check_windowed_matrix( transpose_test_matrix(FULL),                  100 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(throws_if_window_too_thin) {
	constexpr size_t MIN_DIM_VALUE =  2;
	constexpr size_t MAX_VALUE     = 10;
	for (const size_t &smaller_dim_ctr : irange( MIN_DIM_VALUE, MAX_VALUE ) ) {
		for (const size_t &dim_diff_ctr : indices( MAX_VALUE ) ) {
			BOOST_CHECK_THROW(get_window_start_a_for_b__offset_1( smaller_dim_ctr,                smaller_dim_ctr + dim_diff_ctr, dim_diff_ctr, 1 ), invalid_argument_exception);
			BOOST_CHECK_THROW(get_window_start_a_for_b__offset_1( smaller_dim_ctr + dim_diff_ctr, smaller_dim_ctr,                dim_diff_ctr, 1 ), invalid_argument_exception);
			BOOST_CHECK_THROW(get_window_stop_a_for_b__offset_1(  smaller_dim_ctr,                smaller_dim_ctr + dim_diff_ctr, dim_diff_ctr, 1 ), invalid_argument_exception);
			BOOST_CHECK_THROW(get_window_stop_a_for_b__offset_1(  smaller_dim_ctr + dim_diff_ctr, smaller_dim_ctr ,               dim_diff_ctr, 1 ), invalid_argument_exception);
		}
	}
}

/// \brief This testcase was useful for uncovering that replicating previous behaviour meant sometimes producing negative matrix indices!
BOOST_AUTO_TEST_CASE(troublesome_test_case) {
	BOOST_CHECK_EQUAL( 16_z, get_window_upper_and_lower_part_widths(     59, 20, 71        ).first  );
	BOOST_CHECK_EQUAL( 54_z, get_window_upper_and_lower_part_widths(     59, 20, 71        ).second );
	BOOST_CHECK_EQUAL( -4,   get_window_indexing_offset_for_b__offset_1( 59, 20, 71,    20 )        );
	BOOST_CHECK_EQUAL( -3,   get_window_matrix_a_index__offset_1(        59, 20, 71, 1, 20 )        );
	BOOST_CHECK_EQUAL(  4_z, get_window_start_a_for_b__offset_1(         59, 20, 71,    20 )        );
	BOOST_CHECK_EQUAL( 59_z, get_window_stop_a_for_b__offset_1(          59, 20, 71,    20 )        );
}

/// \brief Check that check_indices_are_within_window() behaves as expected on some previously seen data
BOOST_AUTO_TEST_CASE(check_indices_are_within_window_works) {
	BOOST_CHECK_NO_THROW_DIAG( check_indices_are_within_window(9, 7, 11, 8, 2) );
}

BOOST_AUTO_TEST_SUITE_END()

