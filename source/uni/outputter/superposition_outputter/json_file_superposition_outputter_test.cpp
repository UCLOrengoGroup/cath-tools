/// \file
/// \brief The json_file_superposition_outputter test suite

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

#include <boost/filesystem/path.hpp>
#include <boost/test/unit_test.hpp>

#include "chopping/region/region.hpp"
#include "common/file/temp_file.hpp"
#include "outputter/superposition_outputter/json_file_superposition_outputter.hpp"
#include "test/predicate/files_equal.hpp"
#include "test/superposition_fixture.hpp"

using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

using boost::filesystem::path;

namespace cath {
	namespace test {

		/// \brief Assist in testing json_file_superposition_outputter
		struct json_file_superposition_outputter_test_suite_fixture : protected superposition_fixture {
		protected:
			~json_file_superposition_outputter_test_suite_fixture() noexcept = default;

			const temp_file temp_json_file{ "cath_tools_test_temp_file.json_file_superposition_outputter.%%%%" };
			const path temp_json_filename = get_filename( temp_json_file );

			ostringstream err_ss;

			const string json_file_base = "1c0pA01_1hdoA00";
			const string suffix         = ".sup_json";

			const path expected_pretty_file  = TEST_SUP_JSON_DIR() / ( json_file_base + ".pretty"  + suffix );
			const path expected_compact_file = TEST_SUP_JSON_DIR() / ( json_file_base + ".compact" + suffix );

		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(json_file_superposition_outputter_test_suite, cath::test::json_file_superposition_outputter_test_suite_fixture)

BOOST_AUTO_TEST_CASE(writes_correct_pretty_style_json_file_without_errors) {
	json_file_superposition_outputter{ temp_json_filename, json_style::PRETTY }.output_superposition( the_sup_con, err_ss, "dummy_name" );

	BOOST_CHECK_FILES_EQUAL             ( temp_json_filename, expected_pretty_file );
	BOOST_CHECK( err_ss.str().empty() );
}


BOOST_AUTO_TEST_CASE(writes_correct_compact_style_json_file_without_errors) {
	json_file_superposition_outputter{ temp_json_filename, json_style::COMPACT }.output_superposition( the_sup_con, err_ss, "dummy_name" );

	BOOST_CHECK_FILES_EQUAL             ( temp_json_filename, expected_compact_file );
	BOOST_CHECK( err_ss.str().empty() );
}


BOOST_AUTO_TEST_SUITE_END()

