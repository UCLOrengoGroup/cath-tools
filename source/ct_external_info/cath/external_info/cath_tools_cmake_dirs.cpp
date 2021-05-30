/// \file
/// \brief The cath_tools_cmake_dirs definitions

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

#include <string_view>

#include "cath/external_info/cath_tools_cmake_dirs.hpp"
#include "cath/external_info/cath_tools_cmake_dirs_impl.hpp"

using ::std::string_view;

/// The CMAKE_BINARY_DIR (ie base build dir) of the CMake run that is building this
string_view cath::cath_tools_cmake_binary_dir() {
	return CMAKE_BINARY_DIR;
}

/// The CMAKE_SOURCE_DIR (ie base source dir) of the CMake run that is building this
string_view cath::cath_tools_cmake_source_dir() {
	return CMAKE_SOURCE_DIR;
}
