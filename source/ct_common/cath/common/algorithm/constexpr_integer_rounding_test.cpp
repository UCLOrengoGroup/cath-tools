/// \file
/// \brief The integer rounding test suite

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

#include <boost/test/unit_test.hpp>

#include "cath/common/algorithm/constexpr_integer_rounding.hpp"
#include "cath/common/size_t_literal.hpp"

using namespace cath::common;

/// \brief TODOCUMENT
BOOST_AUTO_TEST_SUITE(constexpr_integer_rounding_test_suite)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	static_assert( round_div_up  (  3_z ,  3_z ) ==  1_z, "" );
	static_assert( round_div_up  (  3   ,  3   ) ==  1  , "" );
	static_assert( round_div_up  (  4_z ,  3_z ) ==  2_z, "" );
	static_assert( round_div_up  (  4   ,  3   ) ==  2  , "" );
	static_assert( round_div_up  (  5_z ,  3_z ) ==  2_z, "" );
	static_assert( round_div_up  (  5   ,  3   ) ==  2  , "" );
	static_assert( round_div_up  (  6_z ,  3_z ) ==  2_z, "" );
	static_assert( round_div_up  (  6   ,  3   ) ==  2  , "" );

	static_assert( round_down_mod(  3_z ,  3_z ) ==  3_z, "" );
	static_assert( round_down_mod(  3   ,  3   ) ==  3  , "" );
	static_assert( round_down_mod(  4_z ,  3_z ) ==  3_z, "" );
	static_assert( round_down_mod(  4   ,  3   ) ==  3  , "" );
	static_assert( round_down_mod(  5_z ,  3_z ) ==  3_z, "" );
	static_assert( round_down_mod(  5   ,  3   ) ==  3  , "" );
	static_assert( round_down_mod(  6_z ,  3_z ) ==  6_z, "" );
	static_assert( round_down_mod(  6   ,  3   ) ==  6  , "" );

	static_assert( round_up_mod  (  3_z ,  3_z ) ==  3_z, "" );
	static_assert( round_up_mod  (  3   ,  3   ) ==  3  , "" );
	static_assert( round_up_mod  (  4_z ,  3_z ) ==  6_z, "" );
	static_assert( round_up_mod  (  4   ,  3   ) ==  6  , "" );
	static_assert( round_up_mod  (  5_z ,  3_z ) ==  6_z, "" );
	static_assert( round_up_mod  (  5   ,  3   ) ==  6  , "" );
	static_assert( round_up_mod  (  6_z ,  3_z ) ==  6_z, "" );
	static_assert( round_up_mod  (  6   ,  3   ) ==  6  , "" );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
