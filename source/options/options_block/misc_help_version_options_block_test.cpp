/// \file
/// \brief The misc_help_version_options_block test suite

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

#include "common/exception/invalid_argument_exception.hpp"
#include "options/options_block/misc_help_version_options_block.hpp"
#include "options/options_block/options_block_tester.hpp"

using namespace cath::opts;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The misc_help_version_options_block_test_suite_fixture to assist in testing misc_help_version_options_block
		struct misc_help_version_options_block_test_suite_fixture : protected options_block_tester {
		protected:
			~misc_help_version_options_block_test_suite_fixture() noexcept = default;

			misc_help_version_options_block the_options_block;
			const string IGNORE_OPT = { "positional-that-should-be-ignored" };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(misc_help_version_options_block_test_suite, cath::test::misc_help_version_options_block_test_suite_fixture)

/// \brief Check that if nothing is requested then that is handled correctly
BOOST_AUTO_TEST_CASE(handles_nothing_requested) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT }
	);
	BOOST_CHECK_EQUAL( false, the_options_block.get_help()    );
	BOOST_CHECK_EQUAL( false, the_options_block.get_version() );
}

/// \brief Check that if help is requested then that is handled correctly
BOOST_AUTO_TEST_CASE(handles_help_request) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT,
		  "--" + misc_help_version_options_block::PO_HELP }
	);
	BOOST_CHECK_EQUAL( true,  the_options_block.get_help()    );
	BOOST_CHECK_EQUAL( false, the_options_block.get_version() );
}

/// \brief Check that if version is requested then that is handled correctly
BOOST_AUTO_TEST_CASE(handles_version_request) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT,
		  "--" + misc_help_version_options_block::PO_VERSION }
	);
	BOOST_CHECK_EQUAL( false, the_options_block.get_help()    );
	BOOST_CHECK_EQUAL( true,  the_options_block.get_version() );
}

/// \brief Check that if both help and version are requested then that is handled correctly
BOOST_AUTO_TEST_CASE(handles_help_and_version_requested) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT,
		  "--" + misc_help_version_options_block::PO_HELP,
		  "--" + misc_help_version_options_block::PO_VERSION }
	);
	BOOST_CHECK_EQUAL( true, the_options_block.get_help()    );
	BOOST_CHECK_EQUAL( true, the_options_block.get_version() );
}

BOOST_AUTO_TEST_SUITE_END()

