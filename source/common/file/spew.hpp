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

#ifndef _CATH_TOOLS_SOURCE_COMMON_FILE_SPEW_H
#define _CATH_TOOLS_SOURCE_COMMON_FILE_SPEW_H

#include <boost/filesystem/path.hpp>

#include "common/file/open_fstream.hpp"

#include <fstream>

namespace cath {
	namespace common {

		/// \brief Write the specified string to the specified file
		///
		/// This is named after Perl's Path::Class::File spew() method
		inline void spew(const boost::filesystem::path &arg_file,  ///< The file to which the string should be written
		                 const std::string             &arg_string ///< The string to write to the file
		                 ) {
			std::ofstream output_stream;
			open_ofstream( output_stream, arg_file );
			output_stream << arg_string;
			output_stream.close();
		}

		/// \brief Write the specified istream to the specified file
		///
		/// This is named after Perl's Path::Class::File spew() method
		inline void spew(const boost::filesystem::path &arg_file,   ///< The file to which the string should be written
		                 std::istream                  &arg_istream ///< The istream to write to the file
		                 ) {
			std::ofstream output_stream;
			open_ofstream( output_stream, arg_file );
			arg_istream.clear();
			arg_istream.seekg( 0 );
			output_stream << arg_istream.rdbuf();
			output_stream.close();
		}

	} // namespace common
} // namespace cath

#endif
