/// \file
/// \brief The less_than_helper test suite

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

#include <boost/function.hpp>
#include <boost/range/irange.hpp>

#include "common/less_than_helper.h"
#include "common/size_t_literal.h"
#include "test/global_test_constants.h"

#include <tuple>
#include <utility>

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::irange;

namespace cath {
	namespace test {

		/// \brief The less_than_helper_test_suite_fixture to assist in testing less_than_helper
		struct less_than_helper_test_suite_fixture : protected global_test_constants {
		protected:
			~less_than_helper_test_suite_fixture() noexcept = default;

			/// \brief TODOCUMENT
			using test_type     = std::tuple<size_t, size_t, size_t>;

			/// \brief TODOCUMENT
			using test_type_vec = std::vector<test_type>;

			/// \brief TODOCUMENT
			bool less_than_using_helper(const test_type &arg_helper_test_type_a, ///< TODOCUMENT
			                            const test_type &arg_helper_test_type_b  ///< TODOCUMENT
			                            ) const;

			/// \brief TODOCUMENT
			test_type_vec make_examples() const;

		};

		/// \brief TODOCUMENT
		bool less_than_helper_test_suite_fixture::less_than_using_helper(const test_type &arg_helper_test_type_a, ///< TODOCUMENT
		                                                                 const test_type &arg_helper_test_type_b  ///< TODOCUMENT
		                                                                 ) const {
			auto the_helper = make_less_than_helper( arg_helper_test_type_a, arg_helper_test_type_b );
			the_helper.register_comparison_field( [] (const test_type &x) { return get<0>( x ); } );
			the_helper.register_comparison_field( [] (const test_type &x) { return get<1>( x ); } );
			the_helper.register_comparison_field( [] (const test_type &x) { return get<2>( x ); } );
			return final_less_than_result( the_helper );
		}


		/// \brief TODOCUMENT
		less_than_helper_test_suite_fixture::test_type_vec less_than_helper_test_suite_fixture::make_examples() const {
			test_type_vec examples;
			for (const size_t ctr_1 : irange( 4_z, 6_z) ) {
				for (const size_t ctr_2 : irange( 4_z, 6_z) ) {
					for (const size_t ctr_3 : irange( 4_z, 6_z) ) {
						examples.emplace_back( ctr_1, ctr_2, ctr_3 );
					}
				}
			}
			return examples;
		}

	}
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(less_than_helper_test_suite, cath::test::less_than_helper_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	const auto examples = make_examples();
	for (const test_type &entry_a : examples) {
		for (const test_type &entry_b : examples) {
			const bool expected_less_than_result = ( entry_a < entry_b );
			BOOST_CHECK_EQUAL( less_than_using_helper( entry_a, entry_b ), expected_less_than_result );
		}
	}


}

BOOST_AUTO_TEST_SUITE_END()

