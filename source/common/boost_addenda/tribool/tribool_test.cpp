/// \file
/// \brief The tribool's test suite

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

#include "common/boost_addenda/tribool/tribool.hpp"

namespace cath { namespace test { } }

using namespace cath::common;
using namespace cath::test;

using boost::indeterminate;
using boost::tribool;

namespace cath {
	namespace test {

		/// \brief The tribool test_suite_fixture to assist in testing tribools
		struct tribool_test_suite_fixture {
		protected:
			~tribool_test_suite_fixture() noexcept = default;

			/// \brief An example true tribool
			static constexpr tribool eg_tb_true  = true;

			/// \brief An example false tribool
			static constexpr tribool eg_tb_false = false;

			/// \brief An example indeterminate tribool
			static constexpr tribool eg_tb_indet = indeterminate;
		};

	} // namespace test
} // namespace cath

constexpr tribool tribool_test_suite_fixture::eg_tb_true;
constexpr tribool tribool_test_suite_fixture::eg_tb_false;
constexpr tribool tribool_test_suite_fixture::eg_tb_indet;

BOOST_FIXTURE_TEST_SUITE(tribool_test_suite, tribool_test_suite_fixture)

/// \TODO Come Boost >= 1.58.0 (even on Travis-CI), make the tribool functions constexpr and change these back to static_asserts
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK  (   is_true              ( eg_tb_true  ) ); ///< \TODO: Come GCC with a fix for https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80485, turn these back into static_assert()s
	BOOST_CHECK  ( ! is_false             ( eg_tb_true  ) );
	BOOST_CHECK  ( ! is_not_true          ( eg_tb_true  ) ); ///< \TODO: Come GCC with a fix for https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80485, turn these back into static_assert()s
	BOOST_CHECK  (   is_not_false         ( eg_tb_true  ) );
	BOOST_CHECK  ( ! indeterminate        ( eg_tb_true  ) );
	BOOST_CHECK  (   is_not_indeterminate ( eg_tb_true  ) );

	BOOST_CHECK  ( ! is_true              ( eg_tb_false ) );
	BOOST_CHECK  (   is_false             ( eg_tb_false ) ); ///< \TODO: Come GCC with a fix for https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80485, turn these back into static_assert()s
	BOOST_CHECK  (   is_not_true          ( eg_tb_false ) );
	BOOST_CHECK  ( ! is_not_false         ( eg_tb_false ) ); ///< \TODO: Come GCC with a fix for https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80485, turn these back into static_assert()s
	BOOST_CHECK  ( ! indeterminate        ( eg_tb_false ) );
	BOOST_CHECK  (   is_not_indeterminate ( eg_tb_false ) );

	BOOST_CHECK  ( ! is_true              ( eg_tb_indet ) );
	BOOST_CHECK  ( ! is_false             ( eg_tb_indet ) );
	BOOST_CHECK  (   is_not_true          ( eg_tb_indet ) );
	BOOST_CHECK  (   is_not_false         ( eg_tb_indet ) );
	BOOST_CHECK  (   indeterminate        ( eg_tb_indet ) );
	BOOST_CHECK  ( ! is_not_indeterminate ( eg_tb_indet ) );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
