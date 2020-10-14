/// \file
/// \brief The display_colour_list_test test suite

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

#include "cath/display_colour/display_colour.hpp"
#include "cath/display_colour/display_colour_list.hpp"

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::test;

namespace cath {
	namespace test {

		/// \brief The display_colour_list_test_suite_fixture to assist in testing display_colour_list_test
		struct display_colour_list_test_suite_fixture {
		protected:
			~display_colour_list_test_suite_fixture() noexcept = default;

		public:
			static constexpr size_t NUM_COLOURS_IN_DEFAULT_COLOURS_STRING = 23;
		};

	}  // namespace test
}  // namespace cath

constexpr size_t display_colour_list_test_suite_fixture::NUM_COLOURS_IN_DEFAULT_COLOURS_STRING;

BOOST_FIXTURE_TEST_SUITE(display_colour_list_test_suite, display_colour_list_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(default_colour_string_produces_list_with_sensible_properties) {
	const display_colour_list list = make_display_colour_list_from_string(display_colour_list::DEFAULT_COLOURS_STRING);
	BOOST_CHECK_EQUAL( NUM_COLOURS_IN_DEFAULT_COLOURS_STRING, list.size() );

	BOOST_CHECK_EQUAL( list.colour_of_index(  0 ), colour_of_mod_index( list,  0 ) );
	BOOST_CHECK_EQUAL( list.colour_of_index(  1 ), colour_of_mod_index( list,  1 ) );
	BOOST_CHECK_EQUAL( list.colour_of_index( 21 ), colour_of_mod_index( list, 21 ) );
	BOOST_CHECK_EQUAL( list.colour_of_index( 22 ), colour_of_mod_index( list, 22 ) );
	BOOST_CHECK_EQUAL( list.colour_of_index(  0 ), colour_of_mod_index( list, 23 ) );
	BOOST_CHECK_EQUAL( list.colour_of_index(  1 ), colour_of_mod_index( list, 24 ) );

//	BOOST_CHECK_EQUAL( "cath_tools_defined_colour_0",  name_of_colour_of_mod_index(list,  0) );
//	BOOST_CHECK_EQUAL( "cath_tools_defined_colour_1",  name_of_colour_of_mod_index(list,  1) );
//	BOOST_CHECK_EQUAL( "cath_tools_defined_colour_21", name_of_colour_of_mod_index(list, 21) );
//	BOOST_CHECK_EQUAL( "cath_tools_defined_colour_22", name_of_colour_of_mod_index(list, 22) );
//	BOOST_CHECK_EQUAL( "cath_tools_defined_colour_0",  name_of_colour_of_mod_index(list, 23) );
//	BOOST_CHECK_EQUAL( "cath_tools_defined_colour_1",  name_of_colour_of_mod_index(list, 24) );
}

BOOST_AUTO_TEST_SUITE_END()

