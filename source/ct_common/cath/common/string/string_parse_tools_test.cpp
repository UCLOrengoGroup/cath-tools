/// \file
/// \brief The string_parse_tools test suite

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

#include "cath/common/size_t_literal.hpp"
#include "cath/common/string/string_parse_tools.hpp"

#include "cath/test/global_test_constants.hpp"

namespace cath { namespace test { } }

using namespace ::cath::common;
using namespace ::cath::test;

using ::std::string;

namespace cath {
	namespace test {

		/// \brief The string_parse_tool_test_suite_fixture to assist in testing string_parse_tools
		struct string_parse_tool_test_suite_fixture : protected global_test_constants {
		protected:
			~string_parse_tool_test_suite_fixture() noexcept = default;

			const string pdb_line = "ATOM    189  CZ2 TRP A 584       5.401  40.241  -4.793  1.00 10.59";
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(string_parse_tool_test_suite, string_parse_tool_test_suite_fixture)

BOOST_AUTO_TEST_CASE(rejects_invalid_strings) {
	BOOST_CHECK_THROW( parse_double_from_substring( "a12345b",  0, 6 ), invalid_argument_exception );
	BOOST_CHECK_THROW( parse_int_from_substring   ( "a12345b",  0, 6 ), invalid_argument_exception );
	BOOST_CHECK_THROW( parse_ulong_from_substring ( "a12345b",  0, 6 ), invalid_argument_exception );
	BOOST_CHECK_THROW( parse_double_from_substring( "a12345b",  1, 6 ), invalid_argument_exception );
	BOOST_CHECK_THROW( parse_int_from_substring   ( "a12345b",  1, 6 ), invalid_argument_exception );
	BOOST_CHECK_THROW( parse_ulong_from_substring ( "a12345b",  1, 6 ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(rejects_middle_space) {
	BOOST_CHECK_THROW( parse_double_from_substring( " 12 45 ",  0, 7 ), invalid_argument_exception );
	BOOST_CHECK_THROW( parse_int_from_substring   ( " 12 45 ",  0, 7 ), invalid_argument_exception );
	BOOST_CHECK_THROW( parse_ulong_from_substring ( " 12 45 ",  0, 7 ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(accepts_side_space) {
	BOOST_CHECK_EQUAL( parse_double_from_substring( " 12345 ",  0, 7 ), 12345 );
	BOOST_CHECK_EQUAL( parse_int_from_substring   ( " 12345 ",  0, 7 ), 12345 );
	BOOST_CHECK_EQUAL( parse_ulong_from_substring ( " 12345 ",  0, 7 ), 12345 );
}

BOOST_AUTO_TEST_CASE(accepts_correct_length) {
	BOOST_CHECK_EQUAL( parse_double_from_substring(  "12345",   0,  5 ), 12345 );
	BOOST_CHECK_EQUAL( parse_int_from_substring   (  "12345",   0,  5 ), 12345 );
	BOOST_CHECK_EQUAL( parse_ulong_from_substring (  "12345",   0,  5 ), 12345 );
	BOOST_CHECK_EQUAL( parse_double_from_substring( "a12345b",  1,  5 ), 12345 );
	BOOST_CHECK_EQUAL( parse_int_from_substring   ( "a12345b",  1,  5 ), 12345 );
	BOOST_CHECK_EQUAL( parse_ulong_from_substring ( "a12345b",  1,  5 ), 12345 );
}

BOOST_AUTO_TEST_CASE(parses_double) {
	BOOST_CHECK_EQUAL( parse_double_from_substring( pdb_line,  30,  8 ),   5.401 );
	BOOST_CHECK_EQUAL( parse_double_from_substring( pdb_line,  38,  8 ),  40.241 );
	BOOST_CHECK_EQUAL( parse_double_from_substring( pdb_line,  46,  8 ),  -4.793 );
	BOOST_CHECK_EQUAL( parse_double_from_substring( pdb_line,  54,  6 ),   1.00  );
	BOOST_CHECK_EQUAL( parse_double_from_substring( pdb_line,  60,  6 ),  10.59  );
}

BOOST_AUTO_TEST_CASE(parses_int) {
	BOOST_CHECK_EQUAL( parse_int_from_substring   ( pdb_line,  22,  4 ), 584     );
}

BOOST_AUTO_TEST_CASE(parses_unsigned_long_int) {
	BOOST_CHECK_EQUAL( parse_ulong_from_substring ( pdb_line,   6,  5 ), 189_z   );
}

BOOST_AUTO_TEST_CASE(dumb_trim_string_ref_works) {
	const string source = " billy bob  ";
	BOOST_CHECK_EQUAL( dumb_trim_string_ref( source ), "billy bob" );
}

BOOST_AUTO_TEST_SUITE_END()
