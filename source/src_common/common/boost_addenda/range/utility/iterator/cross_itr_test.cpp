/// \file
/// \brief The cross_itr's test suite

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

#include <boost/algorithm/string/trim.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "common/type_aliases.hpp"
#include "test/boost_addenda/boost_check_equal_ranges.hpp"

#include <algorithm>
#include <tuple>

using namespace cath;
using namespace cath::common;
using namespace std;

BOOST_TEST_DONT_PRINT_LOG_VALUE( bool_size_str_tpl     )
BOOST_TEST_DONT_PRINT_LOG_VALUE( bool_size_str_tpl_vec )


namespace std {
	/// \brief Temporary hacky solution to just get the test working for now
	///
	/// It'd make more sense to use BOOST_TEST_DONT_PRINT_LOG_VALUE() but
	/// that doesn't seem to work with BOOST_CHECK_EQUAL_COLLECTIONS (or BOOST_CHECK_EQUAL_RANGES)
	ostream & operator<<(ostream                 &prm_os,   ///< TODOCUMENT
	                     const bool_size_str_tpl &prm_tuple ///< TODOCUMENT
	                     ) {
		ostringstream temp_ss;
		prm_os << "tuple<bool,size_t,string>(";
		prm_os << boolalpha << get<0>( prm_tuple ) << noboolalpha;
		prm_os << ", ";
		prm_os << get<1>( prm_tuple );
		prm_os << ", ";
		prm_os << get<2>( prm_tuple );
		prm_os << ")";
		return prm_os;
	}
}  // namespace std

namespace cath {
	namespace test {

		/// \brief The cross_itr test_suite_fixture to assist in testing cross_itrs
		struct cross_itr_test_suite_fixture {
		protected:
			~cross_itr_test_suite_fixture() noexcept = default;

			/// \brief TODOCUMENT
			const bool_deq false_true        = { false, true, true };

			/// \brief TODOCUMENT
			const size_vec numbers           = { 12, 11, 10 };

			/// \brief TODOCUMENT
			const str_vec  strings           = { "mary.........", "susan.........", "sally........." };

			/// \brief TODOCUMENT
			const bool_deq premod_false_true = { true, false, false };

			/// \brief TODOCUMENT
			const size_vec premod_numbers    = { 3, 2, 1 };

			/// \brief TODOCUMENT
			const str_vec  premod_strings    = { "mary", "susan", "sally" };

			/// \brief TODOCUMENT
			bool_deq mutable_false_true      = premod_false_true;

			/// \brief TODOCUMENT
			size_vec mutable_numbers         = premod_numbers;

			/// \brief TODOCUMENT
			str_vec  mutable_strings         = premod_strings;

			/// \brief TODOCUMENT
			const bool_size_str_tpl_vec crossed = {
				make_tuple( false, 12, "mary........."  ),
				make_tuple( false, 12, "susan........." ),
				make_tuple( false, 12, "sally........." ),
				make_tuple( false, 11, "mary........."  ),
				make_tuple( false, 11, "susan........." ),
				make_tuple( false, 11, "sally........." ),
				make_tuple( false, 10, "mary........."  ),
				make_tuple( false, 10, "susan........." ),
				make_tuple( false, 10, "sally........." ),
				make_tuple( true,  12, "mary........."  ),
				make_tuple( true,  12, "susan........." ),
				make_tuple( true,  12, "sally........." ),
				make_tuple( true,  11, "mary........."  ),
				make_tuple( true,  11, "susan........." ),
				make_tuple( true,  11, "sally........." ),
				make_tuple( true,  10, "mary........."  ),
				make_tuple( true,  10, "susan........." ),
				make_tuple( true,  10, "sally........." ),
				make_tuple( true,  12, "mary........."  ),
				make_tuple( true,  12, "susan........." ),
				make_tuple( true,  12, "sally........." ),
				make_tuple( true,  11, "mary........."  ),
				make_tuple( true,  11, "susan........." ),
				make_tuple( true,  11, "sally........." ),
				make_tuple( true,  10, "mary........."  ),
				make_tuple( true,  10, "susan........." ),
				make_tuple( true,  10, "sally........." )
			};
		};

		struct modifier final {
			void operator()(tuple<bool &, size_t &, string &> x
			                ) const {
				get<0>( x ) = ! get<0>( x );
				get<1>( x ) += 1;
				get<2>( x ) += ".";
			}
		};

	}  // namespace test
}  // namespace cath

/// \brief Test suite to perform basic checks of the cross_itr
BOOST_FIXTURE_TEST_SUITE(cross_itr_test_suite, cath::test::cross_itr_test_suite_fixture)



BOOST_AUTO_TEST_SUITE(iterator_ts)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(itr) {
	const auto begin_itr = make_cross_itr    ( false_true, numbers, strings );
	const auto end_itr   = make_end_cross_itr( false_true, numbers, strings );
	BOOST_CHECK_EQUAL_RANGES( copy_build<bool_size_str_tpl_vec>( begin_itr, end_itr ), crossed );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(itr_to_mutable) {
	auto begin_itr = make_cross_itr    ( mutable_false_true, mutable_numbers, mutable_strings );
	auto end_itr   = make_end_cross_itr( mutable_false_true, mutable_numbers, mutable_strings );
	std::for_each( begin_itr, end_itr, cath::test::modifier() );
	BOOST_CHECK_EQUAL_RANGES( false_true, mutable_false_true );
	BOOST_CHECK_EQUAL_RANGES( numbers,    mutable_numbers    );
	BOOST_CHECK_EQUAL_RANGES( strings,    mutable_strings    );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(range_ts)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(range) {
	const auto the_const_range = cross( false_true, numbers, strings );
	BOOST_CHECK_EQUAL_RANGES( copy_build<bool_size_str_tpl_vec>( the_const_range ), crossed );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(mutable_range) {
	auto the_range = cross( mutable_false_true, mutable_numbers, mutable_strings );
	boost::range::for_each( the_range, cath::test::modifier() );
	BOOST_CHECK_EQUAL_RANGES( false_true, mutable_false_true );
	BOOST_CHECK_EQUAL_RANGES( numbers,    mutable_numbers    );
	BOOST_CHECK_EQUAL_RANGES( strings,    mutable_strings    );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(non_mutable_tuple_ts)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(tuple_of_refs) {
	auto       the_tuple = tie  ( false_true, numbers, strings );
	const auto the_cross = cross( the_tuple );
	BOOST_CHECK_EQUAL_RANGES( copy_build<bool_size_str_tpl_vec>( the_cross ), crossed );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(tuple_of_vals) {
	auto       the_tuple = tuple<const bool_deq, const size_vec, const str_vec>( false_true, numbers, strings );
	const auto the_cross = cross( the_tuple );
	BOOST_CHECK_EQUAL_RANGES( copy_build<bool_size_str_tpl_vec>( the_cross ), crossed );		
}

// Demonstrates that this can successfully use a temporary made by tie
BOOST_AUTO_TEST_CASE(const_tuple_of_refs) {
	const auto the_cross = cross( tie( false_true, numbers, strings ) );
	BOOST_CHECK_EQUAL_RANGES( copy_build<bool_size_str_tpl_vec>( the_cross ), crossed );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(const_tuple_of_vals) {
	const auto the_tuple = tuple<const bool_deq, const size_vec, const str_vec>( false_true, numbers, strings );
	const auto the_cross = cross( the_tuple );
	BOOST_CHECK_EQUAL_RANGES( copy_build<bool_size_str_tpl_vec>( the_cross ), crossed );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(mutable_tuple_ts)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(tuple_of_refs_to_mutable) {
	auto the_tuple = tie  ( mutable_false_true, mutable_numbers, mutable_strings );
	auto the_cross = cross( the_tuple );
	boost::range::for_each( the_cross, cath::test::modifier() );
	BOOST_CHECK_EQUAL_RANGES( false_true, mutable_false_true );
	BOOST_CHECK_EQUAL_RANGES( numbers,    mutable_numbers    );
	BOOST_CHECK_EQUAL_RANGES( strings,    mutable_strings    );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(tuple_of_mutable_vals) {
	auto the_tuple = make_tuple( premod_false_true, premod_numbers, premod_strings );
	auto the_cross = cross     ( the_tuple );
	boost::range::for_each( the_cross, cath::test::modifier() );
	BOOST_CHECK_EQUAL_RANGES( false_true, get<0>( the_tuple ) );
	BOOST_CHECK_EQUAL_RANGES( numbers,    get<1>( the_tuple ) );
	BOOST_CHECK_EQUAL_RANGES( strings,    get<2>( the_tuple ) );
}

/// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(const_tuple_of_refs_to_mutable) {
//	const auto the_cross = cross( tie( mutable_false_true, mutable_numbers, mutable_strings ) );
//	boost::range::for_each( the_cross, cath::test::modifier() );
//	BOOST_CHECK_EQUAL_RANGES( false_true, mutable_false_true );
//	BOOST_CHECK_EQUAL_RANGES( numbers,    mutable_numbers    );
//	BOOST_CHECK_EQUAL_RANGES( strings,    mutable_strings    );
//}

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(const_tuple_of_mutable_vals) {
//	const auto the_tuple = make_tuple( premod_false_true, premod_numbers, premod_strings );
//	auto       the_cross = cross     ( the_tuple );
//	boost::range::for_each( the_cross, cath::test::modifier() );
//	BOOST_CHECK_EQUAL_RANGES( false_true, get<0>( the_tuple ) );
//	BOOST_CHECK_EQUAL_RANGES( numbers,    get<1>( the_tuple ) );
//	BOOST_CHECK_EQUAL_RANGES( strings,    get<2>( the_tuple ) );
//}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE_END()
