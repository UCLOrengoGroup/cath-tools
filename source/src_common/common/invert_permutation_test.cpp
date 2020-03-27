/// \file
/// \brief The invert_permutation test suite

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

#include "invert_permutation.hpp"

#include <boost/test/unit_test.hpp>

#include "test/boost_addenda/boost_check_equal_ranges.hpp"
#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::common;

namespace cath {
	namespace test {

		/// \brief The invert_permutation_test_suite_fixture to assist in testing invert_permutation
		struct invert_permutation_test_suite_fixture : protected global_test_constants {
		protected:
			~invert_permutation_test_suite_fixture() noexcept = default;

			/// \brief A simple permutation, forwards
			const size_vec permutation_forward  = { 4, 1, 0, 3, 2 };

			/// \brief The same simple permutation, backwards
			const size_vec permutation_backward = { 2, 1, 4, 3, 0 };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(invert_permutation_test_suite, cath::test::invert_permutation_test_suite_fixture)

/// \brief Check that inverting a simple forward permutation correctly generates the backward permutation
BOOST_AUTO_TEST_CASE(invert_forward) {
	const size_vec inverted_forward = invert_permutation( permutation_forward );
	BOOST_CHECK_EQUAL_RANGES( permutation_backward, inverted_forward );
}

/// \brief Check that inverting a simple backward permutation correctly generates the forward permutation
BOOST_AUTO_TEST_CASE(invert_backward) {
	const size_vec inverted_backward = invert_permutation( permutation_backward );
	BOOST_CHECK_EQUAL_RANGES( permutation_forward, inverted_backward );
}

BOOST_AUTO_TEST_SUITE_END()
