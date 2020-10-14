/// \file
/// \brief The simple_chopping_format test suite

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

#include "cath/chopping/chopping_format/simple_chopping_format.hpp"

using namespace ::cath;
using namespace ::cath::chop;

BOOST_AUTO_TEST_SUITE(simple_chopping_format_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( simple_chopping_format{}.parse_segment( "121-232[K]"    ), make_simple_region( 'K', 121,      232      ) );
	BOOST_CHECK_EQUAL( simple_chopping_format{}.parse_segment( "1(B)-99(C)[S]" ), make_simple_region( 'S',   1, 'B',  99, 'C' ) );
}

BOOST_AUTO_TEST_SUITE_END()
