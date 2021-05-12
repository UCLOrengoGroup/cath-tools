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

#include <filesystem>
#include <fstream>

using namespace ::cath::common;
using namespace ::cath::common::detail;

using ::std::ifstream;
using ::std::ios;
using ::std::ios_base;
using ::std::ofstream;
using ::std::filesystem::path;

/// \brief Open an existing ifstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit
///
/// \param prm_fstream  The fstream to open
/// \param prm_filename The file to open
/// \param prm_mode     The mode with which to open the file
void cath::common::open_ifstream( ifstream &prm_fstream, const path &prm_filename, const ios_base::openmode &prm_mode ) {
	return open_existing_fstream_impl( prm_fstream, prm_filename, prm_mode );
}

/// \brief Open an existing ofstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit or failbit
///
/// \param prm_fstream  The fstream to open
/// \param prm_filename The file to open
/// \param prm_mode     The mode with which to open the file
void cath::common::open_ofstream( ofstream &prm_fstream, const path &prm_filename, const ios_base::openmode &prm_mode ) {
	return open_existing_fstream_impl( prm_fstream, prm_filename, prm_mode );
}

/// \brief Construct, open and return an ifstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit
///
/// \param prm_filename The file to open
/// \param prm_mode     The mode with which to open the file
ifstream cath::common::open_ifstream( const ::std::filesystem::path &prm_filename, const ios_base::openmode &prm_mode ) {
	return open_fstream_impl<ifstream>( prm_filename, prm_mode );
}

/// \brief Construct, open and return an ofstream with a path, throwing runtime_error_exception on error, and leave exceptions set to throw on badbit or failbit
///
/// \param prm_filename The file to open
/// \param prm_mode     The mode with which to open the file
ofstream cath::common::open_ofstream( const ::std::filesystem::path &prm_filename, const ios_base::openmode &prm_mode ) {
	return open_fstream_impl<ofstream>( prm_filename, prm_mode );
}
