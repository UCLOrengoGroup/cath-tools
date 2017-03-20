/// \file
/// \brief The chain_label test suite

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

#include "chain_label.hpp"

#include "test/global_test_constants.hpp"

using namespace cath;

namespace cath {
	namespace test {

		/// \brief The chain_label_test_suite_fixture to assist in testing chain_label
		struct chain_label_test_suite_fixture: protected global_test_constants {
		protected:
			~chain_label_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(chain_label_test_suite, cath::test::chain_label_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	static_assert(     chain_label( 'a' ) == chain_label( 'a' )  , "" );
	static_assert( ! ( chain_label( 'a' ) == chain_label( 'b' ) ), "" );
	static_assert( ! ( chain_label( 'a' ) == chain_label( 'A' ) ), "" );
	BOOST_CHECK_EQUAL( chain_label( 'a' ).to_string(), "a" );
	BOOST_CHECK_EQUAL( chain_label( 'b' ).to_string(), "b" );
}

BOOST_AUTO_TEST_SUITE_END()
