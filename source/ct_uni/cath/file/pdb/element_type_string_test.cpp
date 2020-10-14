/// \file
/// \brief The element_type_string test suite

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

#include "cath/common/cpp14/make_unique.hpp"
#include "cath/file/pdb/element_type_string.hpp"

#include <memory>

using namespace ::cath;
using namespace ::cath::file;

using ::std::string;

BOOST_AUTO_TEST_SUITE(element_type_string_test_suite)

BOOST_AUTO_TEST_CASE(get_coarse_element_type_works) {
	BOOST_CHECK_EQUAL( get_coarse_element_type( "C"   ), coarse_element_type::CARBON       );
	BOOST_CHECK_EQUAL( get_coarse_element_type( "CA"  ), coarse_element_type::CARBON_ALPHA );
	BOOST_CHECK_EQUAL( get_coarse_element_type( "CB"  ), coarse_element_type::CARBON_BETA  );
	BOOST_CHECK_EQUAL( get_coarse_element_type( "N"   ), coarse_element_type::NITROGEN     );
	BOOST_CHECK_EQUAL( get_coarse_element_type( "O"   ), coarse_element_type::OXYGEN       );

	BOOST_CHECK_EQUAL( get_coarse_element_type( ""    ), coarse_element_type::NON_CORE     );
	BOOST_CHECK_EQUAL( get_coarse_element_type( "CAA" ), coarse_element_type::NON_CORE     );
	BOOST_CHECK_EQUAL( get_coarse_element_type( "BR"  ), coarse_element_type::NON_CORE     );
}

BOOST_AUTO_TEST_CASE(get_coarse_element_type_of_element_type_string) {
	BOOST_CHECK_EQUAL( get_coarse_element_type( element_type_string{ { { ' ', 'C', 'A', ' ' } } } ), coarse_element_type::CARBON_ALPHA );
	BOOST_CHECK_EQUAL( get_coarse_element_type( element_type_string{ { { ' ', 'C', 'A', ' ' } } } ), coarse_element_type::CARBON_ALPHA );
	BOOST_CHECK_EQUAL( get_coarse_element_type( element_type_string{ { { ' ', 'C', ' ', ' ' } } } ), coarse_element_type::CARBON       );
	BOOST_CHECK_EQUAL( get_coarse_element_type( element_type_string{ { { ' ', 'C', 'A', ' ' } } } ), coarse_element_type::CARBON_ALPHA );
	BOOST_CHECK_EQUAL( get_coarse_element_type( element_type_string{ { { ' ', 'C', 'B', ' ' } } } ), coarse_element_type::CARBON_BETA  );
	BOOST_CHECK_EQUAL( get_coarse_element_type( element_type_string{ { { ' ', 'N', ' ', ' ' } } } ), coarse_element_type::NITROGEN     );
	BOOST_CHECK_EQUAL( get_coarse_element_type( element_type_string{ { { ' ', 'O', ' ', ' ' } } } ), coarse_element_type::OXYGEN       );

	BOOST_CHECK_EQUAL( get_coarse_element_type( element_type_string{ { { ' ', 'C', 'A', 'A' } } } ), coarse_element_type::NON_CORE     );
	BOOST_CHECK_EQUAL( get_coarse_element_type( element_type_string{ { { ' ', 'B', 'R', ' ' } } } ), coarse_element_type::NON_CORE     );
}

BOOST_AUTO_TEST_SUITE_END()