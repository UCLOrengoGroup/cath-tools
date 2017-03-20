/// \file
/// \brief The constexpr_floor test suite

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

#include "common/algorithm/constexpr_floor.hpp"

using namespace cath::common;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The constexpr_floor_test_suite_fixture to assist in testing constexpr_floor
		struct constexpr_floor_test_suite_fixture {
		protected:
			~constexpr_floor_test_suite_fixture() noexcept = default;

			static_assert( is_same<decltype( constexpr_floor( declval<double>() ) ), double>::value, "constexpr_floor( double ) should return a double" );
			static_assert( is_same<decltype( constexpr_floor( declval<float >() ) ), float >::value, "constexpr_floor( float  ) should return a float" );

			static_assert( constexpr_floor( -2.00 ) == -2.00, "constexpr_floor( -2.00 ) should return -2.00" );
			static_assert( constexpr_floor( -1.75 ) == -2.00, "constexpr_floor( -1.75 ) should return -2.00" );
			static_assert( constexpr_floor( -1.50 ) == -2.00, "constexpr_floor( -1.50 ) should return -2.00" );
			static_assert( constexpr_floor( -1.25 ) == -2.00, "constexpr_floor( -1.25 ) should return -2.00" );
			static_assert( constexpr_floor( -1.00 ) == -1.00, "constexpr_floor( -1.00 ) should return -1.00" );
			static_assert( constexpr_floor( -0.75 ) == -1.00, "constexpr_floor( -0.75 ) should return -1.00" );
			static_assert( constexpr_floor( -0.50 ) == -1.00, "constexpr_floor( -0.50 ) should return -1.00" );
			static_assert( constexpr_floor( -0.25 ) == -1.00, "constexpr_floor( -0.25 ) should return -1.00" );
			static_assert( constexpr_floor(  0.00 ) ==  0.00, "constexpr_floor(  0.00 ) should return  0.00" );
			static_assert( constexpr_floor(  0.25 ) ==  0.00, "constexpr_floor(  0.25 ) should return  0.00" );
			static_assert( constexpr_floor(  0.50 ) ==  0.00, "constexpr_floor(  0.50 ) should return  0.00" );
			static_assert( constexpr_floor(  0.75 ) ==  0.00, "constexpr_floor(  0.75 ) should return  0.00" );
			static_assert( constexpr_floor(  1.00 ) ==  1.00, "constexpr_floor(  1.00 ) should return  1.00" );
			static_assert( constexpr_floor(  1.25 ) ==  1.00, "constexpr_floor(  1.25 ) should return  1.00" );
			static_assert( constexpr_floor(  1.50 ) ==  1.00, "constexpr_floor(  1.50 ) should return  1.00" );
			static_assert( constexpr_floor(  1.75 ) ==  1.00, "constexpr_floor(  1.75 ) should return  1.00" );
			static_assert( constexpr_floor(  2.00 ) ==  2.00, "constexpr_floor(  2.00 ) should return  2.00" );

		};
	}  // namespace test
}  // namespace cath
