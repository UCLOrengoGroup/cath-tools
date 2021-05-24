/// \file
/// \brief The constexpr_find test suite

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

#include "cath/common/algorithm/constexpr_find.hpp"

//#include "cath/test/global_test_constants.hpp"

using namespace ::cath::common;
using namespace ::std;

namespace {

		/// \brief The constexpr_find_test_suite_fixture to assist in testing constexpr_find
		struct constexpr_find_test_suite_fixture {
		protected:
			~constexpr_find_test_suite_fixture() noexcept = default;

//			 constexpr array<int, 0> empty_array{ { } };
//			 constexpr int empty_find_result = constexpr_find( empty_array, 1 ); ///< Should fail to compile with sensible message

			static constexpr array<pair<size_t, size_t>, 3> single_num_to_num {{ { 3u, 13u } }};

			static_assert( get<1>( constexpr_find   ( single_num_to_num,  3u ) ) == 13u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<1>( constexpr_find<0>( single_num_to_num,  3u ) ) == 13u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<0>( constexpr_find<1>( single_num_to_num, 13u ) ) ==  3u, "Failure in static_assert test of constexpr_find()" );

			static constexpr array<pair<size_t, size_t>, 5> num_to_num {{
				{ 3u, 13u },
				{ 4u, 14u },
				{ 5u, 15u },
				{ 6u, 16u },
				{ 7u, 17u }
			}};

			static_assert( get<1>( constexpr_find   ( num_to_num,  3u ) ) == 13u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<1>( constexpr_find   ( num_to_num,  4u ) ) == 14u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<1>( constexpr_find   ( num_to_num,  5u ) ) == 15u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<1>( constexpr_find   ( num_to_num,  6u ) ) == 16u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<1>( constexpr_find   ( num_to_num,  7u ) ) == 17u, "Failure in static_assert test of constexpr_find()" );

			static_assert( get<1>( constexpr_find<0>( num_to_num,  3u ) ) == 13u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<1>( constexpr_find<0>( num_to_num,  4u ) ) == 14u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<1>( constexpr_find<0>( num_to_num,  5u ) ) == 15u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<1>( constexpr_find<0>( num_to_num,  6u ) ) == 16u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<1>( constexpr_find<0>( num_to_num,  7u ) ) == 17u, "Failure in static_assert test of constexpr_find()" );

			static_assert( get<0>( constexpr_find<1>( num_to_num, 13u ) ) ==  3u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<0>( constexpr_find<1>( num_to_num, 14u ) ) ==  4u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<0>( constexpr_find<1>( num_to_num, 15u ) ) ==  5u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<0>( constexpr_find<1>( num_to_num, 16u ) ) ==  6u, "Failure in static_assert test of constexpr_find()" );
			static_assert( get<0>( constexpr_find<1>( num_to_num, 17u ) ) ==  7u, "Failure in static_assert test of constexpr_find()" );

			// static_assert( get<0>( constexpr_find<1>( num_to_num, 18u ) ) ==  7u, "Failure in static_assert test of constexpr_find()" ); ///< Should fail to compile

		};
} // namespace
