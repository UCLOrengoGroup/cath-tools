/// \file


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

#include "common/c++14/make_unique.h"
#include "file/pdb/element_type_string.h"

#include <memory>

using namespace cath;
using namespace cath::file;

using std::move;
using std::string;

BOOST_AUTO_TEST_SUITE(element_type_string_test_suite)

BOOST_AUTO_TEST_CASE(copy_constructs) {
	auto first_ptr = common::make_unique<element_type_string>( string{ "  me_to_the_bridge  " } );
	element_type_string second{ *first_ptr };
	first_ptr.reset();
	BOOST_CHECK_EQUAL( second.get_element_type(), "me_to_the_bridge" );
}

BOOST_AUTO_TEST_CASE(move_constructs) {
	auto first_ptr = common::make_unique<element_type_string>( string{ "  me_to_the_bridge  " } );
	element_type_string second{ move( *first_ptr ) };
	first_ptr.reset();
	BOOST_CHECK_EQUAL( second.get_element_type(), "me_to_the_bridge" );
}

BOOST_AUTO_TEST_CASE(copy_assigns) {
	auto first_ptr = common::make_unique<element_type_string>( string{ "  me_to_the_bridge  " } );
	element_type_string second{ string{ "  so_long_my_fatal_friend  "} };
	second = *first_ptr;
	first_ptr.reset();
	BOOST_CHECK_EQUAL( second.get_element_type(), "me_to_the_bridge" );
}

BOOST_AUTO_TEST_CASE(move_assigns) {
	auto first_ptr = common::make_unique<element_type_string>( string{ "  me_to_the_bridge  " } );
	element_type_string second{ string{ "  so_long_my_fatal_friend  "} };
	second = move( *first_ptr );
	first_ptr.reset();
	BOOST_CHECK_EQUAL( second.get_element_type(), "me_to_the_bridge" );
}

BOOST_AUTO_TEST_SUITE_END()