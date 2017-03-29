/// \file
/// \brief The string_of_rapidjson_write test suite

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

#include "string_of_rapidjson_write.hpp"

#include <boost/test/auto_unit_test.hpp>

using namespace cath::common;

BOOST_AUTO_TEST_SUITE(string_of_rapidjson_write_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	const auto fn = [] (rapidjson_writer<json_style::PRETTY> &the_writer) {
		the_writer.start_object();
		the_writer.write_key( "a" );
		the_writer.write_null();
		the_writer.end_object();
	};
	BOOST_CHECK_EQUAL( string_of_rapidjson_write< json_style::PRETTY >( fn, 0 ), "{\n"             R"(    "a": null)"             "\n}" );
	BOOST_CHECK_EQUAL( string_of_rapidjson_write< json_style::PRETTY >( fn, 1 ), "{\n"         R"(        "a": null)"         "\n    }" );
	BOOST_CHECK_EQUAL( string_of_rapidjson_write< json_style::PRETTY >( fn, 2 ), "{\n"     R"(            "a": null)"     "\n        }" );
	BOOST_CHECK_EQUAL( string_of_rapidjson_write< json_style::PRETTY >( fn, 3 ), "{\n" R"(                "a": null)" "\n            }" );
}

BOOST_AUTO_TEST_SUITE_END()
