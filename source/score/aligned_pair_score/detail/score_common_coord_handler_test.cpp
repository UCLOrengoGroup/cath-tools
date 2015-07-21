/// \file
/// \brief The score_common_coord_handler test suite

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "score_common_coord_handler.h"

#include <boost/test/unit_test.hpp>

#include "alignment/common_atom_selection_policy/common_atom_select_ca_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.h"
#include "common/boost_check_no_throw_diag.h"
#include "common/test_tools.h"

//#include <iostream>

using namespace cath::align;
using namespace cath::common::test;
using namespace cath::score::detail;
using namespace std;

namespace cath {
	namespace test {

		struct score_common_coord_handler_fixture {
		protected:
			~score_common_coord_handler_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(score_common_coord_handler_test_suite, cath::test::score_common_coord_handler_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(default_ctor_does_not_throw) {
	BOOST_CHECK_NO_THROW_DIAG( score_common_coord_handler{} );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(policy_ctor_does_not_throw) {
	BOOST_CHECK_NO_THROW_DIAG( score_common_coord_handler( common_residue_select_all_policy(), common_atom_select_ca_policy() ) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(copy_ctor_does_not_throw) {
	const score_common_coord_handler the_handler{};
	BOOST_CHECK_NO_THROW_DIAG( score_common_coord_handler{ the_handler } );
}

/// \brief Check the short_suffix_string is correct for a default-constructed score_common_coord_handler
BOOST_AUTO_TEST_CASE(short_suffix_string_is_correct) {
	BOOST_CHECK_EQUAL( score_common_coord_handler().short_suffix_string(), "" );
}

/// \brief Check the long_suffix_string is correct for a default-constructed score_common_coord_handler
BOOST_AUTO_TEST_CASE(long_suffix_string_is_correct) {
	BOOST_CHECK_EQUAL( score_common_coord_handler().long_suffix_string(), "" );
}

/// \brief Check the description_brackets_string is correct for a default-constructed score_common_coord_handler
BOOST_AUTO_TEST_CASE(description_brackets_string_is_correct) {
	BOOST_CHECK_EQUAL( score_common_coord_handler().description_brackets_string(), "" );
}

/// \brief Check operator== and operator!= work over a range of different score_common_coord_handlers
BOOST_AUTO_TEST_CASE(equality_and_inequality_operators_work) {
	// Loop over the possible non-equal pairs of score_common_coord_handlers indices
	check_equality_operators_on_diff_vals_range( get_all_score_common_coord_handlers() );
}

BOOST_AUTO_TEST_SUITE_END()
