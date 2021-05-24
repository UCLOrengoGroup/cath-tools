/// \file


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
#include <boost/test/unit_test.hpp>

#include <boost/numeric/conversion/cast.hpp>

#include "cath/common/argc_argv_faker.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::std;

/// \brief A method to test that argc_argv_faker does what it should for a given list of arguments
static void test_args(const str_vec &prm_args ///< The vector of arguments to be tested
                      ) {
	// Grab the number of arguments
	const auto num_args = prm_args.size();

	// If there are no arguments then check that an attempt to construct a argc_argv_faker
	// throws an invalid_argument_exception.
	if (num_args == 0) {
		BOOST_CHECK_THROW(argc_argv_faker my_argc_argv_faker(prm_args), invalid_argument_exception);
		return;
	}

	// Construct a  argc_argv_faker from the prm_args
	argc_argv_faker my_argc_argv_faker(prm_args);

	// Check the argc_argv_faker gives the correct result for get_argc() and get_argv()
	BOOST_CHECK_EQUAL(num_args, my_argc_argv_faker.get_argc());
	char * * args_ptr = my_argc_argv_faker.get_argv();
	for (const size_t &prm_ctr : indices( num_args ) ) {
		BOOST_CHECK_EQUAL( prm_args[ prm_ctr ], args_ptr[ prm_ctr ] );
	}
	BOOST_CHECK_EQUAL(args_ptr[num_args], static_cast<char *>( nullptr ));
}

BOOST_AUTO_TEST_SUITE(argc_argv_faker_test_suite)

/// \brief Check that argc_argv_faker correctly throws with zero arguments
BOOST_AUTO_TEST_CASE(zero_args) {
	test_args( str_vec() );
}

/// \brief Check that argc_argv_faker correctly handles one argument
BOOST_AUTO_TEST_CASE(one_arg) {
	test_args( {
		"first_arg"
	} );
}

/// \brief Check that argc_argv_faker correctly handles two arguments
BOOST_AUTO_TEST_CASE(two_args) {
	test_args( {
		"first_arg",
		"second_arg"
	} );
}

/// \brief Check that argc_argv_faker correctly handles many arguments
BOOST_AUTO_TEST_CASE(many_args) {
	test_args( {
		"first_arg",
		"second_arg",
		"third_arg",
		"fourth_arg",
		"fifth_arg",
		"sixth_arg",
		"seven_arg",
		"eight_arg",
		"nine_arg",
		"ten_arg"
	} );
}

BOOST_AUTO_TEST_SUITE_END()
