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

#include <boost/test/auto_unit_test.hpp>

#include "chopping/chopping_format/sillitoe_chopping_format.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common;

BOOST_AUTO_TEST_SUITE(sillitoe_chopping_format_test_suite)

BOOST_AUTO_TEST_CASE(throws_on_attempt_to_parse_invalid_segment) {
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( ""      ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "A"     ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "AA"    ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "3:A"   ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "3-:A"  ), invalid_argument_exception );
	BOOST_CHECK_THROW( sillitoe_chopping_format{}.parse_segment( "3-C:A" ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(parses_valid_segment_correctly) {
	// 1qdm D[1qdmC01]1-191:C  D[1qdmC02]192-247:C,103S-104S:C,248-318:C  F319-1318:C
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_segment( ":R"       ), make_simple_region( 'R'                   ) );
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_segment( "7-232:K"  ), make_simple_region( 'K', 7,      232      ) );
	BOOST_CHECK_EQUAL( sillitoe_chopping_format{}.parse_segment( "1B-99C:S" ), make_simple_region( 'S', 1, 'B',  99, 'C' ) );
}

BOOST_AUTO_TEST_SUITE_END()
