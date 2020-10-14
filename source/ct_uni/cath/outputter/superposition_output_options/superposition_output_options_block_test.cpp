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

#include <boost/test/unit_test.hpp>

#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/display/options/display_spec.hpp"
#include "cath/options/options_block/options_block_tester.hpp"
#include "cath/outputter/superposition_output_options/superposition_output_options_block.hpp"
#include "cath/outputter/superposition_outputter/superposition_outputter.hpp"
#include "cath/outputter/superposition_outputter/superposition_outputter_list.hpp"

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
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(superposition_output_options_block_test_suite, cath::test::superposition_output_options_block_test_suite_fixture)

BOOST_AUTO_TEST_CASE(unparsed_has_no_json_file) {
	BOOST_CHECK( superposition_output_options_block{}.get_json_file().empty() );
}

BOOST_AUTO_TEST_CASE(parses_option_for_to_json_file) {
	const auto parsed_block = parse_into_options_block_copy(
		superposition_output_options_block{},
		{ "--" + superposition_output_options_block::PO_SUP_TO_JSON_FILE, "the_filename" }
	);
	BOOST_CHECK_EQUAL( parsed_block.get_json_file(), path( "the_filename" ) );
}

BOOST_AUTO_TEST_CASE(option_for_to_json_file) {
	const auto parsed_block = parse_into_options_block_copy(
		superposition_output_options_block{},
		{ "--" + superposition_output_options_block::PO_SUP_TO_JSON_FILE, "the_filename" }
	);
	const auto outputters = parsed_block.get_superposition_outputters( display_spec{ "", false, false, false, false } );
	BOOST_CHECK_EQUAL( outputters.size(), 1 );
}

BOOST_AUTO_TEST_SUITE_END()

