/// \file
/// \brief The constexpr_is_uniq test suite

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

#include "cath/common/algorithm/constexpr_is_uniq.hpp"

#include "cath/test/global_test_constants.hpp"

#include <array>

using namespace cath::common;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The constexpr_is_uniq_test_suite_fixture to assist in testing constexpr_is_uniq
		struct constexpr_is_uniq_test_suite_fixture : protected global_test_constants {
		protected:
			~constexpr_is_uniq_test_suite_fixture() noexcept = default;

			static constexpr array<int, 0> unique_empty           = { {                  } };

			static constexpr array<int, 1> unique_single          = { { 1                } };

			static constexpr array<int, 2> unique_pair            = { { 1, 2             } };
			static constexpr array<int, 2> non_unique_pair        = { { 1, 1             } };

			static constexpr array<int, 6> unique_multi           = { { 1, 2, 3, 4, 5, 6 } };
			static constexpr array<int, 6> start_match_multi      = { { 1, 1, 3, 4, 5, 6 } };
			static constexpr array<int, 6> mid_match_multi        = { { 1, 2, 3, 3, 5, 6 } };
			static constexpr array<int, 6> end_match_multi        = { { 1, 2, 3, 4, 6, 6 } };
			static constexpr array<int, 6> wide_split_match_multi = { { 1, 2, 3, 4, 5, 1 } };
			static constexpr array<int, 6> split_match_multi      = { { 1, 2, 3, 4, 2, 6 } };

			static_assert(   constexpr_is_uniq( unique_empty           ), "Failure in static_assert() test of constexpr_is_uniq()" );

			static_assert(   constexpr_is_uniq( unique_single          ), "Failure in static_assert() test of constexpr_is_uniq()" );

			static_assert(   constexpr_is_uniq( unique_pair            ), "Failure in static_assert() test of constexpr_is_uniq()" );
			static_assert( ! constexpr_is_uniq( non_unique_pair        ), "Failure in static_assert() test of constexpr_is_uniq()" );

			static_assert(   constexpr_is_uniq( unique_multi           ), "Failure in static_assert() test of constexpr_is_uniq()" );
			static_assert( ! constexpr_is_uniq( start_match_multi      ), "Failure in static_assert() test of constexpr_is_uniq()" );
			static_assert( ! constexpr_is_uniq( mid_match_multi        ), "Failure in static_assert() test of constexpr_is_uniq()" );
			static_assert( ! constexpr_is_uniq( end_match_multi        ), "Failure in static_assert() test of constexpr_is_uniq()" );
			static_assert( ! constexpr_is_uniq( wide_split_match_multi ), "Failure in static_assert() test of constexpr_is_uniq()" );
			static_assert( ! constexpr_is_uniq( split_match_multi      ), "Failure in static_assert() test of constexpr_is_uniq()" );
		};

	}  // namespace test
}  // namespace cath
