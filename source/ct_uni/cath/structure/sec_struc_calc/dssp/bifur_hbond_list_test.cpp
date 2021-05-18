/// \file
/// \brief The bifur_hbond_list test suite

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

#include <boost/test/unit_test.hpp>

#include "bifur_hbond_list.hpp"

namespace cath { namespace test { } }

using namespace ::cath::sec;
using namespace ::cath::test;

using ::std::nullopt;

BOOST_TEST_DONT_PRINT_LOG_VALUE( hbond_half_opt_pair )

namespace cath {
	namespace test {

		/// \brief The bifur_hbond_list_test_suite_fixture to assist in testing bifur_hbond_list
		struct bifur_hbond_list_test_suite_fixture {
		protected:
			~bifur_hbond_list_test_suite_fixture() noexcept = default;

			/// \brief Example hbond partner index
			static constexpr hbond_partner_t index    = 32;

			/// \brief Example hbond partner energy
			static constexpr hbond_energy_t  energy_a = -0.65;
			/// \brief Example hbond partner energy
			static constexpr hbond_energy_t  energy_b = -0.75;
			/// \brief Example hbond partner energy
			static constexpr hbond_energy_t  energy_c = -0.85;

			/// \brief Example hbond_half
			static constexpr hbond_half      a{ index, energy_a };
			/// \brief Example hbond_half
			static constexpr hbond_half      b{ index, energy_b };
			/// \brief Example hbond_half
			static constexpr hbond_half      c{ index, energy_c };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(bifur_hbond_list_test_suite, bifur_hbond_list_test_suite_fixture)

BOOST_AUTO_TEST_CASE(hbond_half_is_constexpr) {
	static_assert( b.index  == index,    "" );
	static_assert( b.energy == energy_b, "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(hbond_half_lower_energy_value_is_bondier_than) {
	static_assert(   is_bondier_than( c, b ), "" );
	static_assert( ! is_bondier_than( b, b ), "" );
	static_assert( ! is_bondier_than( a, b ), "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(hbond_half_pair_updates_correctly) {
	BOOST_CHECK_EQUAL( update_half_bond_pair_copy( { c,       b       }, a ), hbond_half_opt_pair( c, b       ) );
	BOOST_CHECK_EQUAL( update_half_bond_pair_copy( { c,       a       }, b ), hbond_half_opt_pair( c, b       ) );
	BOOST_CHECK_EQUAL( update_half_bond_pair_copy( { b,       a       }, c ), hbond_half_opt_pair( c, b       ) );

	BOOST_CHECK_EQUAL( update_half_bond_pair_copy( { b,       nullopt }, a ), hbond_half_opt_pair( b, a       ) );
	BOOST_CHECK_EQUAL( update_half_bond_pair_copy( { a,       nullopt }, b ), hbond_half_opt_pair( b, a       ) );

	BOOST_CHECK_EQUAL( update_half_bond_pair_copy( { nullopt, nullopt }, b ), hbond_half_opt_pair( b, nullopt ) );
}

BOOST_AUTO_TEST_CASE(bifur_hbond_updates_correctly) {
	BOOST_CHECK_EQUAL( bifur_hbond{}.update_for_this_nh( a ).get_bound_pair_for_this_nh(), hbond_half_opt_pair( a,       nullopt ) );
	BOOST_CHECK_EQUAL( bifur_hbond{}.update_for_this_nh( a ).get_bound_pair_for_this_co(), hbond_half_opt_pair( nullopt, nullopt ) );
	BOOST_CHECK_EQUAL( bifur_hbond{}.update_for_this_co( a ).get_bound_pair_for_this_nh(), hbond_half_opt_pair( nullopt, nullopt ) );
	BOOST_CHECK_EQUAL( bifur_hbond{}.update_for_this_co( a ).get_bound_pair_for_this_co(), hbond_half_opt_pair( a,       nullopt ) );
}

BOOST_AUTO_TEST_CASE(bifur_hbond_list_to_strings_correctly) {
	BOOST_CHECK_EQUAL( to_string( bifur_hbond_list{ 2 } ), "bifur_hbond_list[\n\tbifur_hbond[nh_1st:(             ), nh_2nd(             ), co_1st:(             ), co_2nd:(             )]\n\tbifur_hbond[nh_1st:(             ), nh_2nd(             ), co_1st:(             ), co_2nd:(             )]]" );
}

BOOST_AUTO_TEST_SUITE_END()
