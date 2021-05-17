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

#include <boost/test/unit_test.hpp>

#include "cath/common/boost_addenda/tribool/tribool.hpp"

using namespace ::cath::common;

using ::boost::indeterminate;
using ::boost::tribool;

BOOST_AUTO_TEST_SUITE(tribool_test_suite)

BOOST_AUTO_TEST_CASE( basic ) {
	/// \brief An example true tribool
	constexpr tribool eg_tb_true = true;

	/// \brief An example false tribool
	constexpr tribool eg_tb_false = false;

	/// \brief An example indeterminate tribool
	constexpr tribool eg_tb_indet = indeterminate;

	static_assert(   is_true              ( eg_tb_true  ) );
	static_assert( ! is_false             ( eg_tb_true  ) );
	static_assert( ! is_not_true          ( eg_tb_true  ) );
	static_assert(   is_not_false         ( eg_tb_true  ) );
	static_assert( ! indeterminate        ( eg_tb_true  ) );
	static_assert(   is_not_indeterminate ( eg_tb_true  ) );

	static_assert( ! is_true              ( eg_tb_false ) );
	static_assert(   is_false             ( eg_tb_false ) );
	static_assert(   is_not_true          ( eg_tb_false ) );
	static_assert( ! is_not_false         ( eg_tb_false ) );
	static_assert( ! indeterminate        ( eg_tb_false ) );
	static_assert(   is_not_indeterminate ( eg_tb_false ) );

	static_assert( ! is_true              ( eg_tb_indet ) );
	static_assert( ! is_false             ( eg_tb_indet ) );
	static_assert(   is_not_true          ( eg_tb_indet ) );
	static_assert(   is_not_false         ( eg_tb_indet ) );
	static_assert(   indeterminate        ( eg_tb_indet ) );
	static_assert( ! is_not_indeterminate ( eg_tb_indet ) );
}

BOOST_AUTO_TEST_SUITE_END()
