/// \file
/// \brief The simple_file_read_write test suite

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

#include "cath/common/file/simple_file_read_write.hpp"
#include "cath/common/file/temp_file.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"
#include "cath/test/global_test_constants.hpp"
#include "cath/test/predicate/files_equal.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::std;

namespace {

	using test_tuple_t   = std::tuple<double, size_t, bool, int, string>;
	using test_tuple_vec = std::vector<test_tuple_t>;

	/// \brief The simple_file_read_write_test_suite_fixture to assist in testing simple_file_read_write
	struct simple_file_read_write_test_suite_fixture : protected global_test_constants {
	protected:
		~simple_file_read_write_test_suite_fixture() noexcept = default;

		/// \brief A list of example numbers to test simple file reading / writing
		const doub_vec EXAMPLE_DOUBLES = {
			1.20206,
			1.41421,
			1.61803,
			2.71828,
			3.14159
		};

		/// \brief A list of example tuples to test simple file reading / writing
		const test_tuple_vec EXAMPLE_TUPLES = {
			make_tuple( 2.3, 1_z, true, -2, string( "hello" ) ),
			make_tuple( 2.3, 0_z, true, -1, string( "there" ) ),
			make_tuple( 2.3, 2_z, true,  0, string( "how"   ) ),
			make_tuple( 2.3, 0_z, true,  1, string( "are"   ) ),
			make_tuple( 2.3, 3_z, true,  2, string( "you"   ) )
		};
	};

} // namespace

/// \brief Some basic unit tests for the simple read_file()/write_file() functions
BOOST_FIXTURE_TEST_SUITE(simple_file_read_write_test_suite, simple_file_read_write_test_suite_fixture)

BOOST_AUTO_TEST_SUITE(doubles)

/// \brief Test that reading doubles from an existing file gives the expected result
BOOST_AUTO_TEST_CASE(read_doubles) {
	const auto got = read_file<double>( EXAMPLE_DOUBLES_FILENAME() );
	BOOST_CHECK_EQUAL_RANGES( EXAMPLE_DOUBLES, got );
}

/// \brief Test that writing doubles to a new file creates an identical file
BOOST_AUTO_TEST_CASE(write_doubles) {
	const temp_file temp_write_file(".%%%%-%%%%-%%%%-%%%%.simple-write-file-doubles-test.txt");
	write_file( get_filename( temp_write_file ), EXAMPLE_DOUBLES );
	BOOST_CHECK_FILES_EQUAL( get_filename( temp_write_file ), EXAMPLE_DOUBLES_FILENAME() );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(tuples)

/// \brief Test that reading tuples from an existing file gives the expected result
BOOST_AUTO_TEST_CASE(read_tuples) {
	const auto got = read_file<test_tuple_t>( EXAMPLE_TUPLES_FILENAME() );
	BOOST_TEST( EXAMPLE_TUPLES == got );
}

/// \brief Test that writing tuples to a new file creates an identical file
BOOST_AUTO_TEST_CASE(write_tuple) {
	const temp_file temp_write_file(".%%%%-%%%%-%%%%-%%%%.simple-write-file-tuples-test.txt");
	write_file( get_filename( temp_write_file ), EXAMPLE_TUPLES );
	BOOST_CHECK_FILES_EQUAL( get_filename( temp_write_file ), EXAMPLE_TUPLES_FILENAME() );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE_END()
