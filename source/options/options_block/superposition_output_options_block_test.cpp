/// \file
/// \brief The superposition_output_options_block test suite

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

#include "display/display_spec/display_spec.h"
#include "exception/not_implemented_exception.h"
#include "options/options_block/options_block_tester.h"
#include "options/options_block/superposition_output_options_block.h"
#include "options/outputter/superposition_outputter/superposition_outputter.h"
#include "options/outputter/superposition_outputter/superposition_outputter_list.h"

using namespace cath;
using namespace cath::common;
using namespace cath::opts;

using boost::filesystem::path;

namespace cath {
	namespace test {

		/// \brief The superposition_output_options_block_test_suite_fixture to assist in testing superposition_output_options_block
		struct superposition_output_options_block_test_suite_fixture : protected options_block_tester {
		protected:
			~superposition_output_options_block_test_suite_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(superposition_output_options_block_test_suite, cath::test::superposition_output_options_block_test_suite_fixture)

BOOST_AUTO_TEST_CASE(unparsed_has_no_json_file) {
	BOOST_CHECK( superposition_output_options_block{}.get_json_file().empty() );
}

BOOST_AUTO_TEST_CASE(parses_option_for_to_json_file) {
	const auto parsed_block = parse_into_options_block_copy(
		superposition_output_options_block{},
		{ "--sup-to-json-file", "the_filename" }
	);
	BOOST_CHECK_EQUAL( parsed_block.get_json_file(), path( "the_filename" ) );
}

BOOST_AUTO_TEST_CASE(option_for_to_json_file_not_yet_implemented_outputter) {
	const auto parsed_block = parse_into_options_block_copy(
		superposition_output_options_block{},
		{ "--sup-to-json-file", "the_filename" }
	);
	BOOST_REQUIRE_THROW( parsed_block.get_superposition_outputters( display_spec{ "", false, false, false, false } ), not_implemented_exception );
	BOOST_WARN_MESSAGE( false, "Currently testing for not_implemented_exception on attempt to use JSON superposition outputter" );
}

BOOST_AUTO_TEST_SUITE_END()

