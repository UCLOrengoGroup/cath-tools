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

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/adaptor/adjacented.hpp"
#include "cath/common/boost_addenda/range/adaptor/equal_grouped.hpp"
#include "cath/common/boost_addenda/range/adaptor/lexical_casted.hpp"
#include "cath/common/boost_addenda/range/adaptor/limited.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "cath/test/boost_test_print_type.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::std;

using ::boost::adaptors::filtered;
using ::boost::sub_range;

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

namespace {

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

} // namespace

/// \brief Test suite to perform basic checks of the adaptors
BOOST_FIXTURE_TEST_SUITE(adaptor_test_suite, adaptor_test_suite_fixture)

/// \brief Check that lexical_casted correctly converts between string/size_t range values
BOOST_AUTO_TEST_CASE(lexical_casted_works) {
	BOOST_TEST( copy_build<str_vec>( numbers | lexical_casted<string>() ) == strings );
	BOOST_TEST( copy_build<size_vec>( strings | lexical_casted<size_t>() ) == numbers );
}

/// \brief Check that adjacented correctly returns pairs of references to pairs of consecutive values
BOOST_AUTO_TEST_CASE(adjacented_works) {
	BOOST_TEST( copy_build<size_size_pair_vec>( numbers | adjacented ) == adjacented_numbers );
}

BOOST_AUTO_TEST_CASE(adjacented_of_irange_works) {
	BOOST_TEST( copy_build<size_size_pair_vec>( indices( 4_z ) | adjacented )
	            == ( size_size_pair_vec{ { 0, 1 }, { 1, 2 }, { 2, 3 } } ) );
}

BOOST_AUTO_TEST_CASE(adjacented_of_filtered_irange_works) {
	BOOST_TEST( copy_build<size_size_pair_vec>( indices( 5_z ) | filtered( []( const int &x ) { return ( x % 2 == 0 ); } ) | adjacented )
	            == ( size_size_pair_vec{ { 0, 2 }, { 2, 4 } } ) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(equal_grouped_works) {
	const auto size_vec_of_size_vec_sub_range = []( const sub_range<const size_vec> &x ) {
		return size_vec{ cbegin( x ), cend( x ) };
	};

	BOOST_TEST( transform_build<size_vec_vec>( equivalents | equal_grouped(), size_vec_of_size_vec_sub_range )
	            == grouped_equivalents );
	BOOST_TEST( transform_build<size_vec_vec>( sorted_equivalents | equal_grouped( less<>() ), size_vec_of_size_vec_sub_range )
	            == grouped_sorted_equivalents );
	BOOST_TEST( transform_build<size_vec_vec>( equivalents | equal_grouped(), size_vec_of_size_vec_sub_range )
	            == grouped_equivalents );
	BOOST_TEST( transform_build<size_vec_vec>( sorted_equivalents | equal_grouped( less<>() ), size_vec_of_size_vec_sub_range )
	            == grouped_sorted_equivalents );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(equal_grouped_works_with_sorted) {
	const auto size_vec_of_size_vec_sub_range = []( const sub_range<const size_vec> &x ) {
		return size_vec{ cbegin( x ), cend( x ) };
	};

	BOOST_TEST( transform_build<size_vec_vec>( sorted_equivalents | equal_grouped( less<>() ), size_vec_of_size_vec_sub_range )
	            == grouped_sorted_equivalents );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(limited_works_0) {
	BOOST_TEST( copy_build<size_vec>( numbers | limited( 0_z ) ) == size_vec() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(limited_works_4) {
	BOOST_TEST( copy_build<size_vec>( numbers | limited( 4_z ) ) == first_four_numbers );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(limited_works_size_or_above) {
	BOOST_TEST( copy_build<size_vec>( numbers | limited( numbers.size() ) ) == numbers );
	BOOST_TEST( copy_build<size_vec>( numbers | limited( numbers.size() + 1 ) ) == numbers );
	BOOST_TEST( copy_build<size_vec>( numbers | limited( numbers.size() + 30 ) ) == numbers );
}

BOOST_AUTO_TEST_SUITE_END()
