/// \file
/// \brief The open_fstream class definitions

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

#include "open_fstream.h"

#include <fstream>

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::common::detail;
using namespace std;

/// \brief Open an ifstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit
///
/// \todo Come C++11 and GCC v5.0, make this a proper ifstream factory function
///       (not currently possible because ifstream should be movable in C++11 but this
///        isn't implemented in GCC v4.9)
void cath::common::open_ifstream(ifstream                 &arg_ifstream, ///< TODOCUMENT
                                 const path               &arg_filename, ///< TODOCUMENT
                                 const ios_base::openmode &arg_mode      ///< TODOCUMENT
                                 ) {
	if ( ! exists( arg_filename ) ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			"Cannot open file \"" + arg_filename.string() + "\" for reading because it doesn't exist"
		));
	}

	open_fstream_impl( arg_ifstream, arg_filename, arg_mode, fstream_type::READING );
	arg_ifstream.exceptions( ios::badbit );
}

/// \brief Open an ofstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit or failbit
///
/// \todo Come C++11 and GCC v5.0, make this a proper ofstream factory function
///       (not currently possible because ofstream should be movable in C++11 but this
///        isn't implemented in GCC v4.9)
void cath::common::open_ofstream(ofstream                 &arg_ofstream, ///< TODOCUMENT
                                 const path               &arg_filename, ///< TODOCUMENT
                                 const ios_base::openmode &arg_mode      ///< TODOCUMENT
                                 ) {
	return open_fstream_impl( arg_ofstream, arg_filename, arg_mode, fstream_type::WRITING );
}
