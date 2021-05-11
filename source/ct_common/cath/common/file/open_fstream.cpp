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

#include "open_fstream.hpp"

#include <fstream>

using namespace ::cath::common;
using namespace ::cath::common::detail;
using namespace ::std;

using ::boost::filesystem::path;

/// \brief Open an ifstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit
///
/// \todo Come C++11 and GCC v5.0, make this a proper ifstream factory function
///       (not currently possible because ifstream should be movable in C++11 but this
///        isn't implemented in GCC v4.9)
void cath::common::open_ifstream(ifstream                 &prm_ifstream, ///< TODOCUMENT
                                 const path               &prm_filename, ///< TODOCUMENT
                                 const ios_base::openmode &prm_mode      ///< TODOCUMENT
                                 ) {
	if ( ! exists( prm_filename ) ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			"Cannot open file \"" + prm_filename.string() + "\" for reading because it doesn't exist"
		));
	}

	open_fstream_impl( prm_ifstream, prm_filename, prm_mode, fstream_type::READING );
	prm_ifstream.exceptions( ios::badbit );
}

/// \brief Open an ofstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit or failbit
///
/// \todo Come C++11 and GCC v5.0, make this a proper ofstream factory function
///       (not currently possible because ofstream should be movable in C++11 but this
///        isn't implemented in GCC v4.9)
void cath::common::open_ofstream(ofstream                 &prm_ofstream, ///< TODOCUMENT
                                 const path               &prm_filename, ///< TODOCUMENT
                                 const ios_base::openmode &prm_mode      ///< TODOCUMENT
                                 ) {
	open_fstream_impl( prm_ofstream, prm_filename, prm_mode, fstream_type::WRITING );
}

/// \brief Open an ifstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit
///
/// \param prm_filename  TODOCUMENT
/// \param prm_mode      TODOCUMENT
ifstream cath::common::open_ifstream( const ::std::filesystem::path &prm_filename, const ios_base::openmode &prm_mode ) {
	return open_fstream_impl<ifstream>( prm_filename, prm_mode );
}

/// \brief Open an ofstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit or failbit
///
/// \param prm_filename  TODOCUMENT
/// \param prm_mode      TODOCUMENT
ofstream cath::common::open_ofstream( const ::std::filesystem::path &prm_filename, const ios_base::openmode &prm_mode ) {
	return open_fstream_impl<ofstream>( prm_filename, prm_mode );
}
