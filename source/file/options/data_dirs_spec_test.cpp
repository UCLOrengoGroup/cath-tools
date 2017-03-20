/// \file
/// \brief The data_dirs_spec test suite

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

#include <boost/filesystem/path.hpp>

#include "common/boost_addenda/test/boost_check_equal_ranges.hpp"
#include "common/file/find_file.hpp"
#include "common/file/open_fstream.hpp"
#include "common/file/temp_file.hpp"
#include "common/type_aliases.hpp"
#include "file/options/data_dirs_spec.hpp"

#include <fstream>

using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The data_dirs_spec_test_suite_fixture to assist in testing data_dirs_spec
		struct data_dirs_spec_test_suite_fixture {
		protected:
			~data_dirs_spec_test_suite_fixture() noexcept = default;

			const string   EXAMPLE_PATH_STR    = ":/usr::/var:";
			const path_vec EXAMPLE_DIRECTORIES = { "/usr", "/var" };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(data_dirs_spec_test_suite, cath::test::data_dirs_spec_test_suite_fixture)

/// \brief Check that split_path_into_directories() works as expected on a simple list
///        of directories (with a few spurious colons thrown in for fun)
BOOST_AUTO_TEST_CASE(split_path_into_directories_works) {
	const path_vec got_paths = split_path_into_directories(EXAMPLE_PATH_STR);
	BOOST_CHECK_EQUAL_RANGES( EXAMPLE_DIRECTORIES, got_paths );
}

/// \brief Check that find_file() correctly returns an empty path if it
///        cannot find the requested file in the specified list of directories
BOOST_AUTO_TEST_CASE(find_file_works) {
	// Create a temporary file with some text in it
	const temp_file temp_file("cath_tools_test_temp_file.data_dirs_spec.%%%%");
	const path      temp_file_filename = get_filename( temp_file );
	ofstream temp_file_stream;
	open_ofstream(temp_file_stream, temp_file_filename);
	temp_file_stream << "some text\n";
	temp_file_stream.close();

	const string full_path = temp_file_filename.parent_path().string() + ":" + EXAMPLE_PATH_STR;
	const string basename  = path(temp_file_filename.filename()).string();

	BOOST_CHECK_EQUAL(
		temp_file_filename,
		find_file(
			split_path_into_directories( full_path ),
			basename
		)
	);
}

/// \brief Check that find_file() correctly returns an empty path if it
///        cannot find the requested file in the specified list of directories
BOOST_AUTO_TEST_CASE(find_file_returns_empty_when_nothing_found) {
	BOOST_CHECK_EQUAL(
		path(),
		find_file(EXAMPLE_DIRECTORIES, "file_that_should_be_impossible_to_file")
	);
}

/// \todo Write tests to check:
///        - data_dirs_spec itself
///        - get_path_of_data_file() and find_file() that uses a data_dirs_spec
///        - that find_file() correctly prefers earlier directories (if the requested file's in several)

BOOST_AUTO_TEST_SUITE_END()

