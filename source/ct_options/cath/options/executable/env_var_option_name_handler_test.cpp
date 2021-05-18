/// \file
/// \brief The env_var_option_name_handler test suite

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

#include "cath/options/executable/env_var_option_name_handler.hpp"

using namespace ::cath::opts;
using namespace ::std;

using ::boost::program_options::options_description;

namespace cath {
	namespace test {

		/// \brief The env_var_option_name_handler_test_suite_fixture to assist in testing env_var_option_name_handler
		struct env_var_option_name_handler_test_suite_fixture {
		  protected:
			~env_var_option_name_handler_test_suite_fixture() noexcept = default;

		  public:
			[[nodiscard]] options_description get_options_description() const {
				options_description temp_od;
				temp_od.add_options()( RECOGNISED_OPTION.c_str(), "a dummy options description" );
				return temp_od;
			}

			const string PREFIX               = { "MY_PREFIX_"                   };
			const string CATH_BIN_PREFIX      = { "CATH_TOOLS_"               };
			const string RECOGNISED_ENV_VAR   = { PREFIX + "RECOGNISED_OPTION"   };
			const string RECOGNISED_OPTION    = { "recognised-option"            };
			const string UNRECOGNISED_ENV_VAR = { PREFIX + "UNRECOGNISED_OPTION" };
			const string UNRECOGNISED_OPTION  = { "unrecognised-option"          };

			const options_description         OPTS               = get_options_description();
			const env_var_option_name_handler TEST_HANDLER_TRUE  = { PREFIX, true,  OPTS };
			const env_var_option_name_handler TEST_HANDLER_FALSE = { PREFIX, false, OPTS };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(env_var_option_name_handler_test_suite, cath::test::env_var_option_name_handler_test_suite_fixture)

/// \brief Check that the operator() always returns empty strings when the prefix is unrecognised
BOOST_AUTO_TEST_CASE(returns_empty_with_mismatching_prefix) {
	BOOST_CHECK_EQUAL("", TEST_HANDLER_TRUE ("DOES_NOT_BEGIN_WITH_PREFIX") );
	BOOST_CHECK_EQUAL("", TEST_HANDLER_FALSE("DOES_NOT_BEGIN_WITH_PREFIX") );
}

/// \brief Check that the operator() always returns recognised options that have the correct prefix
BOOST_AUTO_TEST_CASE(allows_recognised_with_matching_prefix) {
	BOOST_CHECK_EQUAL(RECOGNISED_OPTION, TEST_HANDLER_TRUE ( RECOGNISED_ENV_VAR ) );
	BOOST_CHECK_EQUAL(RECOGNISED_OPTION, TEST_HANDLER_FALSE( RECOGNISED_ENV_VAR ) );
}

/// \brief Check that the operator() returns unrecognised options that have the correct prefix
///        only if allow_unknown was set to false (so that parse_environment() will see it and complain)
BOOST_AUTO_TEST_CASE(allows_or_not_unrecognised_with_matching_prefix) {
	BOOST_CHECK_EQUAL("",                  TEST_HANDLER_TRUE ( UNRECOGNISED_ENV_VAR ) );
	BOOST_CHECK_EQUAL(UNRECOGNISED_OPTION, TEST_HANDLER_FALSE( UNRECOGNISED_ENV_VAR ) );
}

/// \brief Check environment_variable_prefix_of_program_name() on a simple example
BOOST_AUTO_TEST_CASE(environment_variable_prefix_of_program_name_works) {
	BOOST_CHECK_EQUAL("TEST_PROGRAM_NAME_", environment_variable_prefix_of_program_name("tesT-progRam-name"));
}

/// \brief Check option_of_environment_variable_and_prefix() on a simple example
BOOST_AUTO_TEST_CASE(option_of_environment_variable_and_prefix_works) {
	BOOST_CHECK_EQUAL("this-is-bob", option_of_environment_variable_and_prefix(PREFIX+"THIS_IS_BOB",        PREFIX         ));
	BOOST_CHECK_EQUAL("dssp-path",   option_of_environment_variable_and_prefix(CATH_BIN_PREFIX+"DSSP_PATH", CATH_BIN_PREFIX));
}

/// \brief Check option_of_environment_variable_and_prefix() returns an empty string if the prefix doesn't match
BOOST_AUTO_TEST_CASE(option_of_environment_variable_and_prefix_rejects_wrong_prefix) {
	BOOST_CHECK_EQUAL("", option_of_environment_variable_and_prefix("MY_OTHER_PREFIX_THIS_IS_BOB",     PREFIX         ));
	BOOST_CHECK_EQUAL("", option_of_environment_variable_and_prefix("SOME_OTHER_ENVIRONMENT_VARIABLE", CATH_BIN_PREFIX));
}

BOOST_AUTO_TEST_SUITE_END()

