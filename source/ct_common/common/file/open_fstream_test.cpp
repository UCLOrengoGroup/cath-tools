/// \file
/// \brief The open_fstream test suite

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

#include <boost/test/unit_test.hpp>

#include "common/file/open_fstream.hpp"
#include "common/file/temp_file.hpp"
#include "test/global_test_constants.hpp"

#include <fstream>

using namespace cath::common;
using namespace std;

using boost::filesystem::path;

namespace cath {
	namespace test {

		/// \brief The open_fstream_test_suite_fixture to assist in testing open_fstream
		struct open_fstream_test_suite_fixture : protected global_test_constants {
		protected:
			~open_fstream_test_suite_fixture() noexcept = default;

		public:
			const temp_file TEST_OUTPUT_TEMP_FILE { "test_open_fstream_output_file.%%%%-%%%%-%%%%-%%%%" };

			const path      TEST_OUTPUT_FILE      = { get_filename( TEST_OUTPUT_TEMP_FILE ) };
			const path TEST_INPUT_FILE            = { ALIGNMENT_FILE() };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(open_fstream_test_suite, cath::test::open_fstream_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ifstream_basic) {
	ifstream test_ifstream;
	open_ifstream(test_ifstream, TEST_INPUT_FILE);
	BOOST_CHECK(test_ifstream.is_open());
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ofstream_basic) {
	ofstream test_ofstream;
	open_ofstream(test_ofstream, TEST_OUTPUT_FILE);
	BOOST_CHECK(test_ofstream.is_open());
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ifstream_throws_runtime_error_exception_for_nonexistent_file) {
	ifstream test_ifstream;
	BOOST_CHECK_THROW( open_ifstream(test_ifstream, NONEXISTENT_FILE()), runtime_error_exception );
}

BOOST_AUTO_TEST_SUITE_END()

