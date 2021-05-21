/// \file
/// \brief The exception_is_equivalentTest class

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

#include "cath/common/exception/exception_is_equivalent.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"

using namespace ::cath::common;
using namespace ::cath::test;

using ::std::string;
using ::std::string_view;

constexpr string_view DEFAULT_EXCEPTION_STRING = "defaultExceptionString";
constexpr string_view EQUIVALENT_TO_EXCEPTION_WITH_DEFAULT_STRING = "otherExceptionString";

BOOST_AUTO_TEST_SUITE(exception_is_equivalentTestSuite)

BOOST_AUTO_TEST_CASE(TypeCheck) {
	const out_of_range_exception out_of_range_exception1( string{ DEFAULT_EXCEPTION_STRING } );
	const out_of_range_exception out_of_range_exception2( string{ DEFAULT_EXCEPTION_STRING } );
	exception_is_equivalent<out_of_range_exception> equivalent_to_out_of_range_exception1(out_of_range_exception1);

	const invalid_argument_exception invalid_argument_exception1( string{ DEFAULT_EXCEPTION_STRING } );
	const invalid_argument_exception invalid_argument_exception2( string{ DEFAULT_EXCEPTION_STRING } );
	exception_is_equivalent<invalid_argument_exception> equivalent_to_invalid_argument_exception1(invalid_argument_exception1);

	BOOST_CHECK(equivalent_to_out_of_range_exception1(out_of_range_exception1));
	BOOST_CHECK(equivalent_to_out_of_range_exception1(out_of_range_exception2));
	BOOST_CHECK(!equivalent_to_out_of_range_exception1(invalid_argument_exception1));
	BOOST_CHECK(!equivalent_to_out_of_range_exception1(invalid_argument_exception2));

	BOOST_CHECK(equivalent_to_invalid_argument_exception1(invalid_argument_exception1));
	BOOST_CHECK(equivalent_to_invalid_argument_exception1(invalid_argument_exception2));
	BOOST_CHECK(!equivalent_to_invalid_argument_exception1(out_of_range_exception1));
	BOOST_CHECK(!equivalent_to_invalid_argument_exception1(out_of_range_exception2));
}

BOOST_AUTO_TEST_CASE(StringCheck) {
	const invalid_argument_exception exception_with_default_string( string{ DEFAULT_EXCEPTION_STRING } );
	const invalid_argument_exception exception_with_other_string( string{ EQUIVALENT_TO_EXCEPTION_WITH_DEFAULT_STRING } );
	exception_is_equivalent<invalid_argument_exception> EQUIVALENT_TO_EXCEPTION_WITH_DEFAULT_STRING(exception_with_default_string);
	exception_is_equivalent<invalid_argument_exception> equivalent_to_exception_with_other_string(exception_with_other_string);

	BOOST_CHECK(EQUIVALENT_TO_EXCEPTION_WITH_DEFAULT_STRING(exception_with_default_string));
	BOOST_CHECK(!EQUIVALENT_TO_EXCEPTION_WITH_DEFAULT_STRING(exception_with_other_string));
	BOOST_CHECK(!equivalent_to_exception_with_other_string(exception_with_default_string));
	BOOST_CHECK(equivalent_to_exception_with_other_string(exception_with_other_string));

}

BOOST_AUTO_TEST_SUITE_END()
