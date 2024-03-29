/// \file
/// \brief The slurp header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_SLURP_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_SLURP_HPP

#include <filesystem>
#include <fstream>

#include "cath/common/file/open_fstream.hpp"

namespace cath::common {

	/// \brief Read a string from the specified file
	///
	/// This is named after Perl's Path::Tiny / Path::Class::File slurp() method
	inline std::string slurp(const ::std::filesystem::path &prm_file ///< The file from which the string should be read
	                         ) {
		std::ifstream input_stream = open_ifstream( prm_file );

		std::string result_str;
		input_stream.seekg( 0, std::ios::end );
		result_str.resize( static_cast<size_t>( input_stream.tellg() ) );
		input_stream.seekg( 0, std::ios::beg );
		input_stream.read( &result_str[0], static_cast<ptrdiff_t>( result_str.size() ) );

		input_stream.close();
		return result_str;
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_SLURP_HPP
