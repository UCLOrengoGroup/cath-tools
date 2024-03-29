/// \file
/// \brief The spew header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_SPEW_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_SPEW_HPP

#include <filesystem>
#include <fstream>

#include "cath/common/file/open_fstream.hpp"

namespace cath::common {

	/// \brief Write the specified string to the specified file
	///
	/// This is named after Perl's Path::Class::File spew() method
	inline void spew(const ::std::filesystem::path &prm_file,  ///< The file to which the string should be written
	                 const std::string             &prm_string ///< The string to write to the file
	                 ) {
		std::ofstream output_stream = open_ofstream( prm_file );
		output_stream << prm_string;
		output_stream.close();
	}

	/// \brief Write the specified istream to the specified file
	///
	/// This is named after Perl's Path::Class::File spew() method
	inline void spew(const ::std::filesystem::path &prm_file,   ///< The file to which the string should be written
	                 std::istream                  &prm_istream ///< The istream to write to the file
	                 ) {
		std::ofstream output_stream = open_ofstream( prm_file );
		prm_istream.clear();
		prm_istream.seekg( 0 );
		output_stream << prm_istream.rdbuf();
		output_stream.close();
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_SPEW_HPP
