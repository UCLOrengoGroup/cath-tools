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
#include <boost/test/auto_unit_test.hpp>

#include <boost/numeric/conversion/cast.hpp>

#include "common/argc_argv_faker.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/type_aliases.hpp"
#include "common/exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::numeric_cast;

namespace cath {
	namespace test {

		/// \brief The argc_argv_faker_test_suite_fixture to assist in testing argc_argv_faker
		struct argc_argv_faker_test_suite_fixture {
		protected:
			~argc_argv_faker_test_suite_fixture() noexcept = default;

			void test_args(const str_vec &) const;
		};

	}  // namespace test
}  // namespace cath

/// \brief A method to test that argc_argv_faker does what it should for a given list of arguments
void cath::test::argc_argv_faker_test_suite_fixture::test_args(const str_vec &arg_args ///< The vector of arguments to be tested
                                                               ) const {
	// Grab the number of arguments
	const auto num_args = arg_args.size();

	// If there are no arguments then check that an attempt to construct a argc_argv_faker
	// throws an invalid_argument_exception.
	if (num_args == 0) {
		BOOST_CHECK_THROW(argc_argv_faker my_argc_argv_faker(arg_args), invalid_argument_exception);
		return;
	}

	// Construct a  argc_argv_faker from the arg_args
	argc_argv_faker my_argc_argv_faker(arg_args);

	// Check the argc_argv_faker gives the correct result for get_argc() and get_argv()
	BOOST_CHECK_EQUAL(num_args, my_argc_argv_faker.get_argc());
	char * * args_ptr = my_argc_argv_faker.get_argv();
	for (const size_t &arg_ctr : indices( num_args ) ) {
		BOOST_CHECK_EQUAL( arg_args[ arg_ctr ], args_ptr[ arg_ctr ] );
	}
	BOOST_CHECK_EQUAL(args_ptr[num_args], static_cast<char *>( nullptr ));
}

BOOST_FIXTURE_TEST_SUITE(argc_argv_faker_test_suite, cath::test::argc_argv_faker_test_suite_fixture)

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
