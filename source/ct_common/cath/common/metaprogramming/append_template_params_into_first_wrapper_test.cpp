/// \file
/// \brief The append_template_params_into_first_wrapper test suite

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

#include "append_template_params_into_first_wrapper.hpp"

#include <boost/test/unit_test.hpp>

#include <tuple>
#include <type_traits>

using namespace ::cath::common;

using ::std::is_same;
using ::std::tuple;

BOOST_AUTO_TEST_SUITE(append_template_params_into_first_wrapper_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert(
		is_same<
			append_template_params_into_first_wrapper_t<
				tuple< int, int >,
				tuple< char >,
				tuple< tuple<int> >,
				tuple< double, double >
			>,
			tuple< int, int, char, tuple<int>, double, double >
		>::value,
		"append_template_params_into_first_wrapper_t should append the templates' parameters into one templates list wrapped by the first template's wrapper"
	);
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
