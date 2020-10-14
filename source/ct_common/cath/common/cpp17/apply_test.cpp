/// \file
/// \brief The apply test suite

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

#include "cath/common/cpp17/apply.hpp"

using namespace cath::common;

using std::make_tuple;

namespace cath {
	namespace test {

		static constexpr int fn_a(const int &x, const double &y) { return x + static_cast<int>( y ); }

		/// \brief The apply_test_suite_fixture to assist in testing apply
		struct apply_test_suite_fixture {
		protected:
			~apply_test_suite_fixture() noexcept = default;

			/// \totest type_apply_stepwise

			/// \totest type_construct_and_apply_stepwise

			/// \totest apply_stepwise

			static_assert(
				::cath::common::apply( fn_a, make_tuple( 1, 2.0 ) ) == 3,
				"Failure in static_assert test of tuple apply()"
			);

		};
	}  // namespace test
}  // namespace cath
