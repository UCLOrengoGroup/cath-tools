/// \file
/// \brief The rapidjson_writer test suite

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

#include "rapidjson_writer.hpp"

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <string>

using namespace ::cath::common;
using namespace ::std::literals::string_literals;

using ::std::make_tuple;
using ::std::sqrt;

BOOST_AUTO_TEST_SUITE(rapidjson_writer_test_suite)

BOOST_AUTO_TEST_CASE(simple_compact_works) {
	BOOST_CHECK_EQUAL( rapidjson_writer< json_style::COMPACT >{}.start_object().write_key_value( "lue", 42u ).end_object().get_cpp_string(), R"({"lue":42})" );
}

BOOST_AUTO_TEST_CASE(simple_pretty_works) {
	const auto expected = "{\n" R"(    "lue": 42)" "\n}";
	BOOST_CHECK_EQUAL( rapidjson_writer< json_style::PRETTY >{}.start_object().write_key_value( "lue", 42u ).end_object().get_cpp_string(), expected );
	BOOST_CHECK_EQUAL( rapidjson_writer<                    >{}.start_object().write_key_value( "lue", 42u ).end_object().get_cpp_string(), expected );
}

BOOST_AUTO_TEST_CASE(everything_works) {
	BOOST_CHECK_EQUAL(
		rapidjson_writer< json_style::PRETTY >{}
			.start_object()
				.write_key( "stuff" ).start_array()
					.write_value( sqrt( 2.0 )                   )
					.write_value(                        -123   ) // int32_t
					.write_value( static_cast< int64_t>( -123 ) )
					.write_value( static_cast<uint32_t>(  123 ) )
					.write_value( static_cast<uint64_t>(  123 ) )
					.write_value( "cstring"                     )
					.write_value( "cppstring"s                  )
				.end_array()
				.write_key_value( "aye",  true  )
				.write_key_value( "nay"s, false )
				.write_key( "null" ).write_null()
				.write_key( "raws" ).start_array()
					.write_raw_string ( R"({"lue-a", 42a})" )
					.write_raw_string ( R"({"lue-a", 42b})" )
					.write_raw_string ( R"({"lue-a", 42c})" )
				.end_array()

			.end_object()
			.get_cpp_string(),
		R"({
    "stuff": [
        1.4142135623730952,
        -123,
        -123,
        123,
        123,
        "cstring",
        "cppstring"
    ],
    "aye": true,
    "nay": false,
    "null": null,
    "raws": [
        {"lue-a", 42a},
        {"lue-a", 42b},
        {"lue-a", 42c}
    ]
})"
	);
}


BOOST_AUTO_TEST_SUITE_END()
