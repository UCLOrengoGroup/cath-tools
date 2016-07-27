/// \file
/// \brief The transform_build test suite

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

#include "transform_build.h"

#include <boost/test/unit_test.hpp>

#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "common/type_aliases.h"

using namespace cath;
using namespace cath::common;
using namespace std;

namespace cath {
	namespace test {

		struct transform_build_fixture {
		protected:
			~transform_build_fixture() noexcept = default;

		public:
			const size_deq      zero_to_ten_deque               = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 };
			//const size_size_map zero_to_ten_map                 = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 };
			const size_set      zero_to_ten_set                 = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 };
			const size_vec      zero_to_ten_vec                 = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 };
			const size_deq      zero_to_thirty_step_three_deque = { 0,  3,  6,  9, 12, 15, 18, 21, 24, 27, 30 };
			//const size_size_map zero_to_thirty_step_three_map   = { 0,  3,  6,  9, 12, 15, 18, 21, 24, 27, 30 };
			const size_set      zero_to_thirty_step_three_set   = { 0,  3,  6,  9, 12, 15, 18, 21, 24, 27, 30 };
			const size_vec      zero_to_thirty_step_three_vec   = { 0,  3,  6,  9, 12, 15, 18, 21, 24, 27, 30 };
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(transform_build_test_suite, cath::test::transform_build_fixture)

/// \brief TODOCUMENT
///
/// \todo Check transform_build() works for (at least) set, map, vector and deque
BOOST_AUTO_TEST_CASE(basic) {
	const size_vec got_zero_to_thirty_step_three = transform_build<size_vec>(
		zero_to_ten_vec,
		[] (const size_t &x) { return x * 3; }
	);
	BOOST_CHECK_EQUAL_RANGES( zero_to_thirty_step_three_vec, got_zero_to_thirty_step_three );
//	BOOST_WARN_MESSAGE( false, "Check that binary version of transform_build() cannot overrun rng1 or rng2" );
//	BOOST_WARN_MESSAGE( false, "Check that unary/binary version of transform_build() works to/from set"     );
//	BOOST_WARN_MESSAGE( false, "Check that unary/binary version of transform_build() works to/from map"     );
//	BOOST_WARN_MESSAGE( false, "Check that unary/binary version of transform_build() works to/from vector"  );
//	BOOST_WARN_MESSAGE( false, "Check that unary/binary version of transform_build() works to/from deque"   );
}

BOOST_AUTO_TEST_SUITE_END()
