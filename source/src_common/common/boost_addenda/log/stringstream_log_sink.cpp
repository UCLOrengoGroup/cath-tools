/// \file
/// \brief The stringstream_log_sink class definitions

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

#include "stringstream_log_sink.hpp"

using namespace cath;

using ::std::ostringstream;
using ::std::string;

/// \brief Get (a non-const reference to) the stringstream of grabbed log output
ostringstream & stringstream_log_sink::stringstream() {
	return out_ss;
}

/// \brief Get (a const reference to) the stringstream of grabbed log output
const ostringstream & stringstream_log_sink::stringstream() const {
	return out_ss;
}

/// \brief Grab the string from the stringstream of grabbed log output
string stringstream_log_sink::str() const {
	return out_ss.str();
}

/// \brief Return whether the stream's string is current empty
bool stringstream_log_sink::str_is_empty() const {
	return out_ss.str().empty();
}
