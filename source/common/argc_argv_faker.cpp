/// \file
/// \brief The argc_argv_faker class definitions

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

#include "argc_argv_faker.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/shared_array.hpp>

#include "common/boost_addenda/range/indices.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::numeric_cast;
using boost::shared_array;

/// \brief A common initialisation function for the two ctors to call
///
/// \todo In C++11, one ctor can call another directly so, once our build moves
///       to the C++11 standard, move this code back into the str_vec ctor
///       and have the other ctor call that one.
void argc_argv_faker::init(const str_vec &arg_arguments
                           ) {
	if (arg_arguments.empty()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to construct argc_argv_faker with no arguments"));
	}

	arguments.reserve(arg_arguments.size());
	for (const string &argument : arg_arguments) {
		shared_array<char> array_pointer(new char[argument.length() + 1]);
		strncpy(array_pointer.get(), argument.c_str(), argument.length() + 1);
		arguments.push_back(array_pointer);
	}

	// \todo Use a C++11 lambda here (when the build moves to the C++11 standard)
	argument_ptrs.reserve(arguments.size());
	for (shared_array<char> &argument : arguments) {
		argument_ptrs.push_back(argument.get());
	}
	argc = numeric_cast<int>(argument_ptrs.size());
	argument_ptrs.push_back( nullptr );
}

/// \brief Constructor for argc_argv_faker from str_vec arguments.
argc_argv_faker::argc_argv_faker(const str_vec &arg_arguments
                                 ) {
	init(arg_arguments);
}

/// \brief Constructor for argc_argv_faker from const argv-style arguments.
argc_argv_faker::argc_argv_faker(const int          &arg_argc,
                                 const char * const  arg_argv[]
                                 ) {
	str_vec argument_strings;
	argument_strings.reserve( numeric_cast<size_t>( arg_argc ) );
	for (const int &arg_ctr : indices( arg_argc ) ) {
		argument_strings.push_back( arg_argv[ arg_ctr ] );
	}
	init( argument_strings );
}

/// \brief Get the argc-like of the data with which this argc_argv_faker was constructed
int & argc_argv_faker::get_argc() {
	return argc;
}

/// \brief Get the argc-like of the data with which this argc_argv_faker was constructed
const int & argc_argv_faker::get_argc() const {
	return argc;
}

/// \brief Get the argv-like of the data with which this argc_argv_faker was constructed
char * * argc_argv_faker::get_argv() {
	return &argument_ptrs.front();
}

/// \brief Get the argv-like of the data with which this argc_argv_faker was constructed
char * const * argc_argv_faker::get_argv() const {
	return &argument_ptrs.front();
}

/// \brief Simple insertion operator for argc_argv_faker
///
/// \relates argc_argv_faker
ostream & cath::operator<<(ostream               &arg_os,             ///< The ostream to which to output the argc_argv_faker
                           const argc_argv_faker &arg_argc_argv_faker ///< The argc_argv_faker to output
                           ) {
	const int            &argc = arg_argc_argv_faker.get_argc();
	const char * const *  argv = arg_argc_argv_faker.get_argv();
	for (const size_t &arg_ctr : indices( numeric_cast<size_t>(argc) ) ) {
		arg_os << ((arg_ctr > 0) ? ", " : "");
		arg_os << argv[arg_ctr];
	}

	return arg_os;
}
