/// \file
/// \brief The temp_file test suite

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

#include <filesystem>
#include <fstream>

#include <boost/test/unit_test.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/file/temp_file.hpp"

using namespace ::cath::common;
using namespace ::std;

using ::std::filesystem::path;

namespace cath {
	namespace test {

		/// \brief The temp_file_test_suite_fixture to assist in testing temp_file
		struct temp_file_test_suite_fixture {
		protected:
			~temp_file_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(temp_file_test_suite, cath::test::temp_file_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(construct_empty) {
	const temp_file empty_temp_file("");
	BOOST_CHECK_EQUAL( false,  has_filename( empty_temp_file ) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	path filename;
	{
		const temp_file basic_temp_file("cath_tools_test_temp_file_%%%%");
		BOOST_CHECK_EQUAL( true, has_filename( basic_temp_file ) );
		filename = get_filename( basic_temp_file );
		BOOST_CHECK( !filename.empty() );
		ofstream file_stream = open_ofstream( filename );
		file_stream << "string\n";
		file_stream.close();
		BOOST_CHECK( exists(filename) );
	}
	BOOST_CHECK( !exists(filename) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ctor_from_full_path_throws) {
	BOOST_CHECK_THROW(temp_file("/tmp/this_file_should_not_get_created.%%%%"), invalid_argument_exception);
}

BOOST_AUTO_TEST_SUITE_END()
