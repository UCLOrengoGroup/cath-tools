/// \file
/// \brief The sillitoe_chopping_format test suite

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

#include "cath/chopping/chopping_format/sillitoe_chopping_format.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common;

using ::std::string;
using ::boost::string_ref;

BOOST_AUTO_TEST_SUITE(sillitoe_chopping_format_test_suite)

BOOST_AUTO_TEST_CASE(throws_on_attempt_to_parse_invalid_segment) {
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( ""      ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "A"     ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "AA"    ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "3:A"   ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "3-:A"  ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "3-C:A" ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(parses_valid_segment) {
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_segment( ":R"       ), make_simple_region( 'R'                   ) );
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_segment( "7-232:K"  ), make_simple_region( 'K', 7,      232      ) );
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_segment( "1B-99C:S" ), make_simple_region( 'S', 1, 'B',  99, 'C' ) );
}

BOOST_AUTO_TEST_CASE(throws_on_attempt_to_parse_invalid_domain) {
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_domain( ""     ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_domain( "D"    ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_domain( "D["   ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_domain( "D[A"  ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_domain( "D[A]" ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(parses_single_segment_domain) {
	const domain expected{ { make_simple_region( 'C', 1, 191 ) } };
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_domain( "1-191:C"  ), expected );
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_domain( "D1-191:C" ), expected );
}

BOOST_AUTO_TEST_CASE(parses_multi_segment_domain) {
	const domain expected{ {
		make_simple_region( 'C', 192,      247      ),
		make_simple_region( 'C', 103, 'S', 104, 'S' ),
		make_simple_region( 'C', 248,      318      ),
	} };
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_domain( "192-247:C,103S-104S:C,248-318:C"  ), expected );
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_domain( "D192-247:C,103S-104S:C,248-318:C" ), expected );
}

BOOST_AUTO_TEST_CASE(parses_domain_name) {
	const domain expected{
		{
			make_simple_region( 'C', 192,      247      ),
			make_simple_region( 'C', 103, 'S', 104, 'S' ),
			make_simple_region( 'C', 248,      318      ),
		},
		"1qdmC02"
	};
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_domain( "D[1qdmC02]192-247:C,103S-104S:C,248-318:C" ), expected );
}

BOOST_AUTO_TEST_CASE(parses_whole_chain_region) {
	const domain expected{ { make_simple_region( 'B' ) }, "1o7iB" };
	BOOST_TEST( sillitoe_chopping_format{}.parse_domain( "D[1o7iB]:B" ) == expected );
}

// Example:
// 1qdm D[1qdmC01]1-191:C  D[1qdmC02]192-247:C,103S-104S:C,248-318:C  F319-1318:C

/// \todo Include functionality to specify an "all" domain

template <typename Fmt>
string parse_and_write_segment(const string_ref &prm_orig_str) {
	const Fmt format{};
	return format.write_region( format.parse_segment( prm_orig_str ) );
}

template <typename Fmt>
string parse_and_write_domain(const string &prm_orig_str) {
	const Fmt format{};
	return format.write_domain( format.parse_domain( prm_orig_str ) );
}

BOOST_AUTO_TEST_CASE(parsed_and_writes_back_to_orig) {
	BOOST_TEST( parse_and_write_segment<sillitoe_chopping_format>( ":R"                                        ) ==  ":R"                                       );
	BOOST_TEST( parse_and_write_segment<sillitoe_chopping_format>( "7-232:K"                                   ) ==  "7-232:K"                                  );
	BOOST_TEST( parse_and_write_segment<sillitoe_chopping_format>( "1B-99C:S"                                  ) ==  "1B-99C:S"                                 );

	BOOST_TEST( parse_and_write_domain <sillitoe_chopping_format>( "D1-191:C"                                  ) == "D1-191:C"                                  );
	BOOST_TEST( parse_and_write_domain <sillitoe_chopping_format>( "D192-247:C,103S-104S:C,248-318:C"          ) == "D192-247:C,103S-104S:C,248-318:C"          );
	BOOST_TEST( parse_and_write_domain <sillitoe_chopping_format>( "D1-191:C"                                  ) == "D1-191:C"                                  );
	BOOST_TEST( parse_and_write_domain <sillitoe_chopping_format>( "D192-247:C,103S-104S:C,248-318:C"          ) == "D192-247:C,103S-104S:C,248-318:C"          );
	BOOST_TEST( parse_and_write_domain <sillitoe_chopping_format>( "D[1o7iB]:B"                                ) == "D[1o7iB]:B"                                );
	BOOST_TEST( parse_and_write_domain <sillitoe_chopping_format>( "D[1qdmC02]192-247:C,103S-104S:C,248-318:C" ) == "D[1qdmC02]192-247:C,103S-104S:C,248-318:C" );
}

BOOST_AUTO_TEST_SUITE_END()
