/// \file
/// \brief The constexpr_modulo_fns test suite

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

//#include <boost/log/trivial.hpp> // ***** TEMPORARY *****
#include <boost/range/combine.hpp>
#include <boost/range/irange.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/algorithm/cxx11/none_of.hpp>

#include "cath/common/algorithm/constexpr_modulo_fns.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/common/type_aliases.hpp"

//#include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::algorithm::none_of;
using boost::irange;
using boost::range::combine;

namespace cath {
	namespace test {

		/// \brief The _test_suite_fixture to assist in testing 
		struct constexpr_modulo_fns_test_suite_fixture {
		protected:
			~constexpr_modulo_fns_test_suite_fixture() noexcept = default;

			/// \brief TODOCUMENT
			void check_crcp(const size_t &prm_result_a, ///< TODOCUMENT
			                const size_t &prm_result_b, ///< TODOCUMENT
			                const size_t &prm_index_a,  ///< TODOCUMENT
			                const size_t &prm_index_b,  ///< TODOCUMENT
			                const size_t &prm_mod_a,    ///< TODOCUMENT
			                const size_t &prm_mod_b     ///< TODOCUMENT
			                ) {
				// Check that the result does indeed hit one of the representative pairs
				BOOST_CHECK_EQUAL( prm_result_a % prm_mod_a, 0_z );
				BOOST_CHECK_EQUAL( prm_result_b % prm_mod_b, 0_z );

				// Check that the result is on the same diagonal as the original point
				// and is at least as far along it as the original point
				BOOST_REQUIRE_GE ( prm_result_a, prm_index_a );
				BOOST_REQUIRE_GE ( prm_result_b, prm_index_b );
				BOOST_CHECK_EQUAL( prm_result_a - prm_index_a, prm_result_b - prm_index_b );

				// Check that there isn't any earlier match (from the original point onwards)
				// that hits a representative pair
				const auto try_range_a = irange( prm_index_a, prm_result_a );
				const auto try_range_b = irange( prm_index_b, prm_result_b );
				BOOST_CHECK( none_of(
					combine( try_range_a, try_range_b ),
					[&] (const boost::tuple<size_t, size_t> &x) {
						const auto &try_a = get<0>( x );
						const auto &try_b = get<1>( x );
						const bool correct_result = ( try_a % prm_mod_a == 0 && try_b % prm_mod_b == 0 );
//						cerr << "Checking none_of [" << try_a << ", " << try_b << "] - result is " << boolalpha << correct_result << noboolalpha << "\n";
						return correct_result;
					}
				) );
			}
		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(constexpr_modulo_fns_test_suite, cath::test::constexpr_modulo_fns_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(constexpr_lcm_works) {
	static_assert( constexpr_lcm( 1_z, 1_z ) ==  1, "constexpr_lcm( 1, 1 ) should be equal to  1" );
	static_assert( constexpr_lcm( 1_z, 2_z ) ==  2, "constexpr_lcm( 1, 2 ) should be equal to  2" );
	static_assert( constexpr_lcm( 1_z, 3_z ) ==  3, "constexpr_lcm( 1, 3 ) should be equal to  3" );
	static_assert( constexpr_lcm( 1_z, 4_z ) ==  4, "constexpr_lcm( 1, 4 ) should be equal to  4" );
	static_assert( constexpr_lcm( 1_z, 5_z ) ==  5, "constexpr_lcm( 1, 5 ) should be equal to  5" );
	static_assert( constexpr_lcm( 1_z, 6_z ) ==  6, "constexpr_lcm( 1, 6 ) should be equal to  6" );
	static_assert( constexpr_lcm( 2_z, 1_z ) ==  2, "constexpr_lcm( 2, 1 ) should be equal to  2" );
	static_assert( constexpr_lcm( 2_z, 2_z ) ==  2, "constexpr_lcm( 2, 2 ) should be equal to  2" );
	static_assert( constexpr_lcm( 2_z, 3_z ) ==  6, "constexpr_lcm( 2, 3 ) should be equal to  6" );
	static_assert( constexpr_lcm( 2_z, 4_z ) ==  4, "constexpr_lcm( 2, 4 ) should be equal to  4" );
	static_assert( constexpr_lcm( 2_z, 5_z ) == 10, "constexpr_lcm( 2, 5 ) should be equal to 10" );
	static_assert( constexpr_lcm( 2_z, 6_z ) ==  6, "constexpr_lcm( 2, 6 ) should be equal to  6" );
	static_assert( constexpr_lcm( 3_z, 1_z ) ==  3, "constexpr_lcm( 3, 1 ) should be equal to  3" );
	static_assert( constexpr_lcm( 3_z, 2_z ) ==  6, "constexpr_lcm( 3, 2 ) should be equal to  6" );
	static_assert( constexpr_lcm( 3_z, 3_z ) ==  3, "constexpr_lcm( 3, 3 ) should be equal to  3" );
	static_assert( constexpr_lcm( 3_z, 4_z ) == 12, "constexpr_lcm( 3, 4 ) should be equal to 12" );
	static_assert( constexpr_lcm( 3_z, 5_z ) == 15, "constexpr_lcm( 3, 5 ) should be equal to 15" );
	static_assert( constexpr_lcm( 3_z, 6_z ) ==  6, "constexpr_lcm( 3, 6 ) should be equal to  6" );
	static_assert( constexpr_lcm( 4_z, 1_z ) ==  4, "constexpr_lcm( 4, 1 ) should be equal to  4" );
	static_assert( constexpr_lcm( 4_z, 2_z ) ==  4, "constexpr_lcm( 4, 2 ) should be equal to  4" );
	static_assert( constexpr_lcm( 4_z, 3_z ) == 12, "constexpr_lcm( 4, 3 ) should be equal to 12" );
	static_assert( constexpr_lcm( 4_z, 4_z ) ==  4, "constexpr_lcm( 4, 4 ) should be equal to  4" );
	static_assert( constexpr_lcm( 4_z, 5_z ) == 20, "constexpr_lcm( 4, 5 ) should be equal to 20" );
	static_assert( constexpr_lcm( 4_z, 6_z ) == 12, "constexpr_lcm( 4, 6 ) should be equal to 12" );
	static_assert( constexpr_lcm( 5_z, 1_z ) ==  5, "constexpr_lcm( 5, 1 ) should be equal to  5" );
	static_assert( constexpr_lcm( 5_z, 2_z ) == 10, "constexpr_lcm( 5, 2 ) should be equal to 10" );
	static_assert( constexpr_lcm( 5_z, 3_z ) == 15, "constexpr_lcm( 5, 3 ) should be equal to 15" );
	static_assert( constexpr_lcm( 5_z, 4_z ) == 20, "constexpr_lcm( 5, 4 ) should be equal to 20" );
	static_assert( constexpr_lcm( 5_z, 5_z ) ==  5, "constexpr_lcm( 5, 5 ) should be equal to  5" );
	static_assert( constexpr_lcm( 5_z, 6_z ) == 30, "constexpr_lcm( 5, 6 ) should be equal to 30" );
	static_assert( constexpr_lcm( 6_z, 1_z ) ==  6, "constexpr_lcm( 6, 1 ) should be equal to  6" );
	static_assert( constexpr_lcm( 6_z, 2_z ) ==  6, "constexpr_lcm( 6, 2 ) should be equal to  6" );
	static_assert( constexpr_lcm( 6_z, 3_z ) ==  6, "constexpr_lcm( 6, 3 ) should be equal to  6" );
	static_assert( constexpr_lcm( 6_z, 4_z ) == 12, "constexpr_lcm( 6, 4 ) should be equal to 12" );
	static_assert( constexpr_lcm( 6_z, 5_z ) == 30, "constexpr_lcm( 6, 5 ) should be equal to 30" );
	static_assert( constexpr_lcm( 6_z, 6_z ) ==  6, "constexpr_lcm( 6, 6 ) should be equal to  6" );
	BOOST_CHECK( true );
}


/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(constexpr_gcd_works) {
	static_assert( constexpr_gcd( 0_z, 0_z ) == 0, "constexpr_gcd( 0, 0 ) should be equal to 0" );
	static_assert( constexpr_gcd( 0_z, 1_z ) == 1, "constexpr_gcd( 0, 1 ) should be equal to 1" );
	static_assert( constexpr_gcd( 0_z, 2_z ) == 2, "constexpr_gcd( 0, 2 ) should be equal to 2" );
	static_assert( constexpr_gcd( 0_z, 3_z ) == 3, "constexpr_gcd( 0, 3 ) should be equal to 3" );
	static_assert( constexpr_gcd( 0_z, 4_z ) == 4, "constexpr_gcd( 0, 4 ) should be equal to 4" );
	static_assert( constexpr_gcd( 0_z, 5_z ) == 5, "constexpr_gcd( 0, 5 ) should be equal to 5" );
	static_assert( constexpr_gcd( 1_z, 0_z ) == 1, "constexpr_gcd( 1, 0 ) should be equal to 1" );
	static_assert( constexpr_gcd( 1_z, 1_z ) == 1, "constexpr_gcd( 1, 1 ) should be equal to 1" );
	static_assert( constexpr_gcd( 1_z, 2_z ) == 1, "constexpr_gcd( 1, 2 ) should be equal to 1" );
	static_assert( constexpr_gcd( 1_z, 3_z ) == 1, "constexpr_gcd( 1, 3 ) should be equal to 1" );
	static_assert( constexpr_gcd( 1_z, 4_z ) == 1, "constexpr_gcd( 1, 4 ) should be equal to 1" );
	static_assert( constexpr_gcd( 1_z, 5_z ) == 1, "constexpr_gcd( 1, 5 ) should be equal to 1" );
	static_assert( constexpr_gcd( 2_z, 0_z ) == 2, "constexpr_gcd( 2, 0 ) should be equal to 2" );
	static_assert( constexpr_gcd( 2_z, 1_z ) == 1, "constexpr_gcd( 2, 1 ) should be equal to 1" );
	static_assert( constexpr_gcd( 2_z, 2_z ) == 2, "constexpr_gcd( 2, 2 ) should be equal to 2" );
	static_assert( constexpr_gcd( 2_z, 3_z ) == 1, "constexpr_gcd( 2, 3 ) should be equal to 1" );
	static_assert( constexpr_gcd( 2_z, 4_z ) == 2, "constexpr_gcd( 2, 4 ) should be equal to 2" );
	static_assert( constexpr_gcd( 2_z, 5_z ) == 1, "constexpr_gcd( 2, 5 ) should be equal to 1" );
	static_assert( constexpr_gcd( 3_z, 0_z ) == 3, "constexpr_gcd( 3, 0 ) should be equal to 3" );
	static_assert( constexpr_gcd( 3_z, 1_z ) == 1, "constexpr_gcd( 3, 1 ) should be equal to 1" );
	static_assert( constexpr_gcd( 3_z, 2_z ) == 1, "constexpr_gcd( 3, 2 ) should be equal to 1" );
	static_assert( constexpr_gcd( 3_z, 3_z ) == 3, "constexpr_gcd( 3, 3 ) should be equal to 3" );
	static_assert( constexpr_gcd( 3_z, 4_z ) == 1, "constexpr_gcd( 3, 4 ) should be equal to 1" );
	static_assert( constexpr_gcd( 3_z, 5_z ) == 1, "constexpr_gcd( 3, 5 ) should be equal to 1" );
	static_assert( constexpr_gcd( 4_z, 0_z ) == 4, "constexpr_gcd( 4, 0 ) should be equal to 4" );
	static_assert( constexpr_gcd( 4_z, 1_z ) == 1, "constexpr_gcd( 4, 1 ) should be equal to 1" );
	static_assert( constexpr_gcd( 4_z, 2_z ) == 2, "constexpr_gcd( 4, 2 ) should be equal to 2" );
	static_assert( constexpr_gcd( 4_z, 3_z ) == 1, "constexpr_gcd( 4, 3 ) should be equal to 1" );
	static_assert( constexpr_gcd( 4_z, 4_z ) == 4, "constexpr_gcd( 4, 4 ) should be equal to 4" );
	static_assert( constexpr_gcd( 4_z, 5_z ) == 1, "constexpr_gcd( 4, 5 ) should be equal to 1" );
	static_assert( constexpr_gcd( 5_z, 0_z ) == 5, "constexpr_gcd( 5, 0 ) should be equal to 5" );
	static_assert( constexpr_gcd( 5_z, 1_z ) == 1, "constexpr_gcd( 5, 1 ) should be equal to 1" );
	static_assert( constexpr_gcd( 5_z, 2_z ) == 1, "constexpr_gcd( 5, 2 ) should be equal to 1" );
	static_assert( constexpr_gcd( 5_z, 3_z ) == 1, "constexpr_gcd( 5, 3 ) should be equal to 1" );
	static_assert( constexpr_gcd( 5_z, 4_z ) == 1, "constexpr_gcd( 5, 4 ) should be equal to 1" );
	static_assert( constexpr_gcd( 5_z, 5_z ) == 5, "constexpr_gcd( 5, 5 ) should be equal to 5" );
	BOOST_CHECK( true );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(chinese_remainder_coprime_pair_works) {
	static_assert(                     detail::extended_euclid_algo(             5_z, 3_z ) == diff_diff_pair( -1,   2 ),  "" );
	static_assert(                     detail::extended_euclid_algo(             3_z, 5_z ) == diff_diff_pair(  2,  -1 ),  "" );
	static_assert(            detail::extended_euclid_algo_products(             5_z, 3_z ) == diff_diff_pair( -5,   6 ),  "" );
	static_assert(            detail::extended_euclid_algo_products(             5_z, 3_z ) == diff_diff_pair( -5,   6 ),  "" );
	static_assert( detail::chinese_remainder_coprime_pair_num_above( 10_z,  4_z, 5_z, 3_z ) == 19_z,                       "" );
	static_assert( detail::chinese_remainder_coprime_pair_num_above(  4_z, 10_z, 3_z, 5_z ) == 19_z,                       "" );

	static_assert(         chinese_remainder_coprime_pair          ( 10_z,  4_z, 5_z, 3_z ) == make_pair( 15_z,  9_z ),    "" );
	static_assert(         chinese_remainder_coprime_pair          (  4_z, 10_z, 3_z, 5_z ) == make_pair(  9_z, 15_z ),    "" );

	static_assert(         chinese_remainder_coprime_pair          (  0_z,  0_z, 1_z, 2_z ) == make_pair(  0_z,  0_z ),    "" );
	static_assert(         chinese_remainder_coprime_pair          (  0_z,  0_z, 2_z, 1_z ) == make_pair(  0_z,  0_z ),    "" );

	static_assert(         chinese_remainder_coprime_pair          (  0_z,  1_z, 1_z, 1_z ) == make_pair(  0_z,  1_z ),    "" );
	static_assert(         chinese_remainder_coprime_pair          (  1_z,  0_z, 1_z, 1_z ) == make_pair(  1_z,  0_z ),    "" );

	static_assert(         chinese_remainder_coprime_pair          (  0_z,  2_z, 2_z, 3_z ) == make_pair(  4_z,  6_z ),    "" );
	static_assert(         chinese_remainder_coprime_pair          (  2_z,  0_z, 3_z, 2_z ) == make_pair(  6_z,  4_z ),    "" );

	// This can ranges can be extended higher for more comprehensive testing
	// (have used 1-10 and 0-100 at the time of writing)
	const auto mod_range   = irange( 1_z,  5_z );
	const auto index_range = indices( 20_z );

	for (const auto &mods_tuple : cross( mod_range, mod_range ) ) {
		const auto &mod_a = get<0>( mods_tuple );
		const auto &mod_b = get<1>( mods_tuple );
		if ( constexpr_gcd( mod_a, mod_b ) == 1 ) {
			for (const auto &indices_tuple : cross( index_range, index_range ) ) {
				const auto &index_a = get<0>( indices_tuple );
				const auto &index_b = get<1>( indices_tuple );
				const auto  result = chinese_remainder_coprime_pair( index_a, index_b, mod_a, mod_b );
				check_crcp( result.first, result.second, index_a, index_b, mod_a, mod_b );
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
