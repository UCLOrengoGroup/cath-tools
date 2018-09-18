/// \file
/// \brief The program_exception_wrapper test file

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

#include <boost/algorithm/string/predicate.hpp>
//#include <boost/scoped_array.hpp>
#include <boost/test/output_test_stream.hpp>

#include "common/argc_argv_faker.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/program_exception_wrapper.hpp"

#include <sstream>
#include <string>

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::test_tools::output_test_stream;
using boost::algorithm::contains;

namespace cath {
	namespace test {

		/// \brief A concrete program_exception_wrapper that throws requested exceptions to allow testing of program_exception_wrapper
		class test_program_exception_wrapper final : public program_exception_wrapper {
			string do_get_program_name() const final {
				return "test_program_exception_wrapper";
			}

			/// \brief Throw the sort of argument requested (or none at all)
			void do_run_program(int argc,
			                    char * argv[]
			                    ) final {
				BOOST_REQUIRE_GE(argc, 1);
				const string first_argument_str( argv[ 0 ] );
				if (first_argument_str == "boost") {
					BOOST_THROW_EXCEPTION(invalid_argument_exception("test_program_exception_wrapper-test-boost-exception"));
				}
				else if (first_argument_str == "std") {
					throw std::exception();
				}
				else if (first_argument_str == "other") {
					throw "test_program_exception_wrapper-test-bare_string-exception";
				}
			}
		};

		/// \brief The program_exception_wrapper_test_suite_fixture to assist in testing program_exception_wrapper
		struct program_exception_wrapper_test_suite_fixture {
		protected:
			~program_exception_wrapper_test_suite_fixture() noexcept = default;

		public:
			void test_arg(const string &prm_string,
			              const string &prm_throw_regexp
			              ) {
				output_test_stream output_stream;
				const str_vec arguments(1, prm_string);
				argc_argv_faker my_argc_argv_faker(arguments);

				test_program_exception_wrapper().run_program(
					my_argc_argv_faker.get_argc(),
					my_argc_argv_faker.get_argv(),
					output_stream
				);
				if ( ! prm_throw_regexp.empty() ) {
					BOOST_CHECK( contains( output_stream.str(), prm_throw_regexp ) );
				}
				else {
					BOOST_CHECK(output_stream.is_empty(true));
				}
			}
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(program_exception_wrapper_test_suite, cath::test::program_exception_wrapper_test_suite_fixture)

/// \brief Test that a boost::exception is caught and described correctly
BOOST_AUTO_TEST_CASE(boost_exception) {
	test_arg("boost", "boost::exception");
}

/// \brief Test that a std::exception is caught and described correctly
BOOST_AUTO_TEST_CASE(std_exception) {
	test_arg("std", "std::exception");
}

/// \brief Test that another exception (here bare string) is caught and described correctly
BOOST_AUTO_TEST_CASE(other_exception) {
	test_arg("other", "exception of unrecognised type");
}

/// \brief Test that all is well if no exception is thrown
BOOST_AUTO_TEST_CASE(no_exception) {
	test_arg("none", "");
}

BOOST_AUTO_TEST_SUITE_END()
