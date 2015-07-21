/// \file
/// \brief The env_var_option_name_handler class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "env_var_option_name_handler.h"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>

using namespace boost::algorithm;
using namespace boost::program_options;
using namespace cath::opts;
using namespace std;

using boost::algorithm::replace_all_copy;
using boost::algorithm::starts_with;

/// \brief Ctor for env_var_option_name_handler
env_var_option_name_handler::env_var_option_name_handler(const string              &arg_prefix,        ///< The prefix string to strip off all environment variable names(eg "CATH_BINARIES_"
                                                         const bool                &arg_allow_unknown, ///< Whether to allow unrecognised options (by not passing them back to parse_environment(), which would complain)
                                                         const options_description &arg_options        ///< The options_description containing the options that should be accepted (can be left to default empty value if arg_allow_unknown is false)
                                                         ) : prefix(arg_prefix),
                                                             allow_unknown(arg_allow_unknown),
                                                             the_options(arg_options) {
}

/// \brief Operator() that takes an environment variable name and converts it to an option name that should be processed or an empty string
///
/// This returns an empty string if:
///  - the specified prefix does not appear at the start of the string or
///  - allow_unknown was specified as true and this option is not recognised by the specified options.
///
/// Otherwise, it returns the result of:
///  - stripping of the prefix,
///  - lower-casing the name and
///  - replacing underscores with hyphens.
string env_var_option_name_handler::operator()(const string &arg_environment_variable_name ///< The environment variable name
                                               ) const {
	const string option_string = option_of_environment_variable_and_prefix(arg_environment_variable_name, prefix);
	// If this name doesn't begin with the prefix, it should just be ignored so return an empty string
	// (note: this test is necessary because an empty option_string shouldn't be passed to find_nothrow())
	if (option_string.empty()) {
		return "";
	}

	// If we're allowing unknown options and this one isn't in the options description,
	// then simply return "" so that Boost Program Options doesn't get a chance to complain
	if (allow_unknown && !the_options.find_nothrow(option_string, false)) {
		return "";
	}

	// Otherwise, pass back the adjusted option name
	return option_string;
}

/// \brief Convert an environment variable to an option name if it begins with the correct prefix, return "" otherwise
///
/// If the string begins with the prefix, this returns the result of:
///  - stripping of the prefix,
///  - lower-casing the name and
///  - replacing underscores with hyphens.
///
/// For example:
///  - the environment variable "CATH_BINARIES_DSSP_PATH" and the prefix "CATH_BINARIES_" would produce "dssp-path"
///  - the environment variable "SOME_OTHER_ENVIRONMENT_VARIABLE" and the prefix "CATH_BINARIES_" would produce ""
///
/// \relates env_var_option_name_handler
string cath::opts::option_of_environment_variable_and_prefix(const string &arg_env_variable_name, ///< The environment variable name
                                                             const string &arg_prefix             ///< The prefix with which the environment variable should begin for it to be accepted
                                                             ) {
	// If this environment variable doesn't start with the correct prefix, ignore it
	if (!starts_with(arg_env_variable_name, arg_prefix)) {
		return "";
	}

	const string unprefixed_name       = erase_first_copy(arg_env_variable_name, arg_prefix);
	const string lower_unprefixed_name = to_lower_copy(unprefixed_name);
	return replace_all_copy(lower_unprefixed_name, "_", "-");
}

/// \brief Convert a program name to an environment variable prefix
///
/// This returns the result of:
///  - upper-casing
///  - replacing all hyphens with underscores
///  - appending a trailing underscore
///
/// For example, "tesT-progRam-name" would produce "TEST_PROGRAM_NAME_"
///
/// \relates env_var_option_name_handler
string cath::opts::environment_variable_prefix_of_program_name(const string &arg_program_name ///< The program name to convert
                                                               ) {
	const string upper_program_name             = to_upper_copy(   arg_program_name            );
	const string underscored_upper_program_name = replace_all_copy(upper_program_name, "-", "_");
	return underscored_upper_program_name + "_";
}
