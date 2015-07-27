/// \file
/// \brief The transform_tuple test suite

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

#include "transform_tuple.h"

#include <boost/test/unit_test.hpp>

#include "common/size_t_literal.h"
#include "common/type_aliases.h"

using namespace cath;
using namespace cath::common;
using namespace std;

BOOST_TEST_DONT_PRINT_LOG_VALUE( size_size_tpl )

namespace cath {
	namespace test {

		struct transform_tuple_fixture {
		protected:
			~transform_tuple_fixture() noexcept = default;

			static constexpr size_size_tpl three_two_tuple = make_tuple( 3_z, 2_z );
			static constexpr size_size_tpl nine_six_tuple  = make_tuple( 9_z, 6_z );
		};

		struct tripler final {
			template <typename... Ts>
			constexpr auto operator()(Ts &... arg_num
			                          ) {
				tie( arg_num... ) = make_tuple( ( 3 * arg_num )... );
			}
		};

		struct tripler_copy final {
			template <typename... Ts>
			constexpr auto operator()(const Ts &... arg_num
			                          ) {
				return make_tuple( ( 3 * arg_num )... );
			}
		};

	}
}

constexpr size_size_tpl cath::test::transform_tuple_fixture::three_two_tuple;
constexpr size_size_tpl cath::test::transform_tuple_fixture::nine_six_tuple;

BOOST_FIXTURE_TEST_SUITE(transform_tuple_test_suite, cath::test::transform_tuple_fixture)

/// \brief Check const transform_tuple() works
BOOST_AUTO_TEST_CASE(const_transform_tuple_works) {
	BOOST_CHECK_EQUAL( 0, 0 );
	static_assert( transform_tuple( three_two_tuple, cath::test::tripler_copy() ) == nine_six_tuple, "Tuple tripling didn't work" );
}

/// \brief Check non-consttransform_tuple() works
BOOST_AUTO_TEST_CASE(non_const_transform_tuple_works) {
	auto non_const_copy = three_two_tuple;
	transform_tuple( non_const_copy, cath::test::tripler() );
	BOOST_CHECK_EQUAL( non_const_copy, nine_six_tuple );
}

BOOST_AUTO_TEST_SUITE_END()
