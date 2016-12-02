/// \file
/// \brief The find_file definitions

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

#include "find_file.hpp"

#include <boost/filesystem.hpp>

using boost::filesystem::exists;
using boost::filesystem::path;
using std::string;

/// \brief Search for a particular file basename through a path of directories
///
/// \returns The found file or an empty path object if one could not be found
path cath::common::find_file(const path_vec &arg_path_dirs, ///< Directories through which to search for the file (in descending order of preference)
                             const string   &arg_basename   ///< The basename of the file to be located (eg 1c0pA01.dssp)
                             ) {
	for (const path &dir : arg_path_dirs) {
		const path potential_file = dir / arg_basename;
		if ( exists( potential_file ) ) {
			return potential_file;
		}
	}
	return path();
}
