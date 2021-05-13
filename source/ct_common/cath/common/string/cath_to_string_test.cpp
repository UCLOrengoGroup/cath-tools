/// \file
/// \brief The cath_to_string test suite

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

#include <optional>
#include <tuple>
#include <utility>

#include <boost/test/unit_test.hpp>

#include "cath/common/size_t_literal.hpp"
#include "cath/common/string/cath_to_string.hpp"

using namespace ::cath::common;

using ::std::nullopt;
using ::std::optional;
using ::std::pair;
using ::std::string;
using ::std::tuple;

BOOST_AUTO_TEST_SUITE( cath_to_string_test_suite )

BOOST_AUTO_TEST_CASE( basic ) {
	BOOST_TEST( cath_to_string( nullopt ) == "nullopt" );
	BOOST_TEST( cath_to_string( optional( 3_z ) ) == "std::optional<unsigned long>(3)" );
	BOOST_TEST( cath_to_string( optional<int>() ) == "std::optional<int>(nullopt)" );

	BOOST_TEST( cath_to_string( pair( 3_z, 2.3 ) ) == "std::pair<unsigned long, double>(3, 2.300000)" );

	BOOST_TEST( cath_to_string( false ) == "false" );
	BOOST_TEST( cath_to_string( true ) == "true" );

	BOOST_TEST( cath_to_string( string{ "hello" } ) == string{ "hello" } );

	BOOST_TEST( cath_to_string( tuple( 3_z, 2.3 ) ) == "std::tuple<unsigned long, double>(3, 2.300000)" );
}

BOOST_AUTO_TEST_SUITE_END()
