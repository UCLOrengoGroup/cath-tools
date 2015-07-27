/// \file
/// \brief The selected_pair test suite

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
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

#include "common/boost_check_no_throw_diag.h"
#include "ssap/selected_pair.h"

using namespace cath;

namespace cath {
	namespace test {

		/// \brief The selected_pair_test_suite_fixture to assist in testing selected_pair
		struct selected_pair_test_suite_fixture {
		protected:
			~selected_pair_test_suite_fixture() noexcept = default;

		public:
			const size_t     TEST_INDEX_A      =  3;
			const size_t     TEST_INDEX_B      =  5;
			const score_type TEST_SCORE        = 10;
			const score_type TEST_HIGHER_SCORE = 20;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(selected_pair_test_suite, cath::test::selected_pair_test_suite_fixture)

/// \brief Check that the basic ctor doesn't throw
BOOST_AUTO_TEST_CASE(ctor_does_not_throw) {
	BOOST_CHECK_NO_THROW_DIAG( selected_pair( TEST_INDEX_A, TEST_INDEX_B, TEST_SCORE ) );
}

/// \brief Check that the getters work as expected
BOOST_AUTO_TEST_CASE(getters_work) {
	const selected_pair test_pair( TEST_INDEX_A, TEST_INDEX_B, TEST_SCORE );
	BOOST_CHECK_EQUAL( TEST_INDEX_A, test_pair.get_index_a() );
	BOOST_CHECK_EQUAL( TEST_INDEX_B, test_pair.get_index_b() );
	BOOST_CHECK_EQUAL( TEST_SCORE,   test_pair.get_score() );
}

/// \brief Check that the less-than operator works as expected
BOOST_AUTO_TEST_CASE(less_than_works) {
	const selected_pair test_pair       ( TEST_INDEX_A, TEST_INDEX_B, TEST_SCORE        );
	const selected_pair test_pair_higher( TEST_INDEX_A, TEST_INDEX_B, TEST_HIGHER_SCORE );

	BOOST_CHECK(     test_pair        < test_pair_higher );
	BOOST_CHECK( ! ( test_pair_higher < test_pair        ) );
	BOOST_CHECK( ! ( test_pair        < test_pair        ) );
	BOOST_CHECK( ! ( test_pair_higher < test_pair_higher ) );
}

BOOST_AUTO_TEST_SUITE_END()

