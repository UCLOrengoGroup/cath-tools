/// \file
/// \brief The adaptors test suite

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

#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/adaptor/adjacented.hpp"
#include "common/boost_addenda/range/adaptor/equal_grouped.hpp"
#include "common/boost_addenda/range/adaptor/lexical_casted.hpp"
#include "common/boost_addenda/range/adaptor/limited.hpp"
#include "common/boost_addenda/test/boost_check_equal_ranges.hpp"
#include "common/boost_addenda/test/boost_check_no_throw_diag.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/pair_insertion_operator.hpp"
#include "common/size_t_literal.hpp"
#include "common/type_aliases.hpp"

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::irange;
using boost::range::for_each;
using boost::sub_range;

namespace std {
	/// \brief Temporary hacky solution to just get the test working for now
	///
	/// It'd make more sense to use BOOST_TEST_DONT_PRINT_LOG_VALUE() but
	/// that doesn't seem to work with BOOST_CHECK_EQUAL_COLLECTIONS (or BOOST_CHECK_EQUAL_RANGES)
	template <typename T>
	ostream & operator<<(ostream         &arg_os, ///< TODOCUMENT
	                     const vector<T> &arg_vec ///< TODOCUMENT
	                     ) {
		ostringstream temp_ss;
		arg_os << "[";
		for (const auto &entry : arg_vec) {
			arg_os << " " << entry;
		}
		arg_os << " ]";
		return arg_os;
	}
}

/// \todo Consider adding a cbegin/cend-like interface for each adaptor, (eg "cequal_grouped") that returns
///       a const equal_grouped interface to a non-const range.
///
/// \todo Consider adding a standard, non-| (ie non-pipe) style interface for each adaptor
///
/// \todo Test mutability of non-const adaptors

// * adjacented;
// * lexical_casted;
// * batched;
// * batched_into_n_batches;
// * equal_grouped;
// * limited;
//   subscripted;
//   iteratored;
//   index_iteratored;
//   permuted;

namespace cath {
	namespace test {

		/// \brief The adaptor test_suite_fixture to assist in testing adaptors
		struct adaptor_test_suite_fixture {
		protected:
			~adaptor_test_suite_fixture() noexcept = default;

			/// \brief TODOCUMENT
			const str_vec            strings            = { "5", "3", "7", "1", "6", "8", "0", "9", "2" };

			/// \brief TODOCUMENT
			const size_vec           numbers            = {  5,   3,   7,   1,   6,   8,   0,   9,   2  };

			/// \brief TODOCUMENT
			const size_vec           first_four_numbers = {  5,   3,   7,   1 };

			/// \brief TODOCUMENT
			const size_size_pair_vec adjacented_numbers = {
				{ 5, 3 },
				{ 3, 7 },
				{ 7, 1 },
				{ 1, 6 },
				{ 6, 8 },
				{ 8, 0 },
				{ 0, 9 },
				{ 9, 2 },
			};

			/// \brief TODOCUMENT
			const size_vec equivalents = { 5, 5, 3, 3, 3, 7, 7, 7 };

			/// \brief TODOCUMENT
			const size_vec_vec grouped_equivalents = {
				{ 5, 5 },
				{ 3, 3, 3 },
				{ 7, 7, 7 }
			};

			/// \brief TODOCUMENT
			const size_vec sorted_equivalents = { 3, 3, 3, 5, 5, 7, 7, 7 };

			/// \brief TODOCUMENT
			const size_vec_vec grouped_sorted_equivalents = {
				{ 3, 3, 3 },
				{ 5, 5 },
				{ 7, 7, 7 }
			};
		};

	}
}

/// \brief Test suite to perform basic checks of the adaptors
BOOST_FIXTURE_TEST_SUITE(adaptor_test_suite, cath::test::adaptor_test_suite_fixture)

/// \brief Check that lexical_casted correctly converts between string/size_t range values
BOOST_AUTO_TEST_CASE(lexical_casted_works) {
	BOOST_CHECK_EQUAL_RANGES( copy_build<str_vec> ( numbers | lexical_casted<string>() ), strings );
	BOOST_CHECK_EQUAL_RANGES( copy_build<size_vec>( strings | lexical_casted<size_t>() ), numbers );
}

/// \brief Check that adjacented correctly returns pairs of references to pairs of consecutive values
BOOST_AUTO_TEST_CASE(adjacented_works) {
	BOOST_CHECK_EQUAL_RANGES(
		copy_build<size_size_pair_vec>( numbers | adjacented ),
		adjacented_numbers
	);
}

BOOST_AUTO_TEST_CASE(adjacented_of_irange_works) {
	BOOST_CHECK_EQUAL_RANGES(
		copy_build<size_size_pair_vec>( irange( 0_z, 4_z ) | adjacented ),
		size_size_pair_vec{ { 0, 1 }, { 1, 2 }, { 2, 3 } }
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(equal_grouped_works) {
	const auto size_vec_of_size_vec_sub_range = [] (const sub_range<const size_vec> &x) {
		return size_vec{ cath::common::cbegin( x ), cath::common::cend( x ) };
	};

	BOOST_CHECK_EQUAL_RANGES(
		transform_build<size_vec_vec>( equivalents        | equal_grouped(                ), size_vec_of_size_vec_sub_range ),
		grouped_equivalents
	);
	BOOST_CHECK_EQUAL_RANGES(
		transform_build<size_vec_vec>( sorted_equivalents | equal_grouped( less<size_t>() ), size_vec_of_size_vec_sub_range ),
		grouped_sorted_equivalents
	);
	BOOST_CHECK_EQUAL_RANGES(
		transform_build<size_vec_vec>( equivalents        | equal_grouped(                ), size_vec_of_size_vec_sub_range ),
		grouped_equivalents
	);
	BOOST_CHECK_EQUAL_RANGES(
		transform_build<size_vec_vec>( sorted_equivalents | equal_grouped( less<size_t>() ), size_vec_of_size_vec_sub_range ),
		grouped_sorted_equivalents
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(equal_grouped_works_with_sorted) {
	const auto size_vec_of_size_vec_sub_range = [] (const sub_range<const size_vec> &x) {
		return size_vec{ cath::common::cbegin( x ), cath::common::cend( x ) };
	};

	BOOST_CHECK_EQUAL_RANGES(
		transform_build<size_vec_vec>( sorted_equivalents | equal_grouped( less<size_t>()  ), size_vec_of_size_vec_sub_range ),
		grouped_sorted_equivalents
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(limited_works_0) {
	BOOST_CHECK_EQUAL_RANGES(
		copy_build<size_vec>( numbers | limited( 0_z ) ),
		size_vec()
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(limited_works_4) {
	BOOST_CHECK_EQUAL_RANGES(
		copy_build<size_vec>( numbers | limited( 4_z ) ),
		first_four_numbers
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(limited_works_size_or_above) {
	BOOST_CHECK_EQUAL_RANGES(
		copy_build<size_vec>( numbers | limited( numbers.size() ) ),
		numbers
	);
	BOOST_CHECK_EQUAL_RANGES(
		copy_build<size_vec>( numbers | limited( numbers.size() + 1 ) ),
		numbers
	);
	BOOST_CHECK_EQUAL_RANGES(
		copy_build<size_vec>( numbers | limited( numbers.size() + 30 ) ),
		numbers
	);
}

BOOST_AUTO_TEST_SUITE_END()
