/// \file
/// \brief The read_string_from_file header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_READ_STRING_FROM_FILE_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_READ_STRING_FROM_FILE_HPP

#include <filesystem>
#include <fstream>

#include "cath/common/file/open_fstream.hpp"

namespace cath::common {

	/// \brief Read the contents of the specified file into a string
	inline std::string read_string_from_file(const ::std::filesystem::path &prm_file ///< The file from which the string should be read
	                                         ) {
		std::ifstream input_stream = open_ifstream( prm_file );

		const std::string result{
			std::istreambuf_iterator<char>( input_stream ),
			std::istreambuf_iterator<char>(              )
		};
		input_stream.close();
		return result;
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_READ_STRING_FROM_FILE_HPP
