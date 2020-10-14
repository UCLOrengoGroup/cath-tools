/// \file
/// \brief The detail_help_options_block test suite

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

#include <boost/test/unit_test.hpp>

#include "common/exception/invalid_argument_exception.hpp"
#include "options/options_block/detail_help_options_block.hpp"
#include "options/options_block/options_block_tester.hpp"

using namespace cath::common;
using namespace cath::opts;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The detail_help_options_block_test_suite_fixture to assist in testing detail_help_options_block
		struct detail_help_options_block_test_suite_fixture : protected options_block_tester {
		protected:
			~detail_help_options_block_test_suite_fixture() noexcept = default;

		public:
			detail_help_options_block the_options_block{ TEST_DESC_AND_HELP_OF_OPTION_NAME() };
			const string              IGNORE_OPT       { "positional-that-should-be-ignored" };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(detail_help_options_block_test_suite, cath::test::detail_help_options_block_test_suite_fixture)

/// \brief Check that if no help is requested then no help is given (and an attempt to retrieve it throws)
BOOST_AUTO_TEST_CASE(handles_no_help_requested) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT }
	);
	BOOST_CHECK_EQUAL( false, the_options_block.has_help_string());
	BOOST_CHECK_THROW( the_options_block.help_string(), invalid_argument_exception );
}

/// \brief Check that if help 1 is requested, then help 1 is given
BOOST_AUTO_TEST_CASE(handles_help_1_requested) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT, "--" + TEST_OPTION_1 }
	);
	BOOST_REQUIRE_EQUAL( true,        the_options_block.has_help_string() );
	BOOST_CHECK_EQUAL(   TEST_HELP_1, the_options_block.help_string()     );
}

/// \brief Check that if help 2 is requested, then help 2 is given
BOOST_AUTO_TEST_CASE(handles_help_2_requested) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT, "--" + TEST_OPTION_2 }
	);
	BOOST_REQUIRE_EQUAL( true,        the_options_block.has_help_string() );
	BOOST_CHECK_EQUAL(   TEST_HELP_2, the_options_block.help_string()     );
}

/// \brief Check that if help 1 and then help 2 is requested, then help 1 is given
BOOST_AUTO_TEST_CASE(handles_help_1_and_2_requested) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT, "--" + TEST_OPTION_1, "--" + TEST_OPTION_2 }
	);
	BOOST_REQUIRE_EQUAL( true,        the_options_block.has_help_string() );
	BOOST_CHECK_EQUAL(   TEST_HELP_1, the_options_block.help_string()     );
}

/// \brief Check that if help 2 and then help 1 is requested, then help 1 is given
BOOST_AUTO_TEST_CASE(handles_help_2_and_1_requested) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT, "--" + TEST_OPTION_2, "--" + TEST_OPTION_1 }
	);
	BOOST_REQUIRE_EQUAL( true,        the_options_block.has_help_string() );
	BOOST_CHECK_EQUAL(   TEST_HELP_1, the_options_block.help_string()     );
}

BOOST_AUTO_TEST_SUITE_END()

