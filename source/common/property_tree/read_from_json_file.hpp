/// \file
/// \brief The read_from_json_file header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_PROPERTY_TREE_READ_FROM_JSON_FILE_H
#define _CATH_TOOLS_SOURCE_COMMON_PROPERTY_TREE_READ_FROM_JSON_FILE_H

#include <boost/filesystem/path.hpp>
// #include <boost/property_tree/json_parser.hpp>

#include "common/file/open_fstream.hpp"
#include "common/property_tree/from_json_string.hpp"

#include <fstream>
#include <string>

namespace cath {
	namespace common {

		/// \brief Build a T from a JSON string (via a ptree)
		///
		/// Requires that there is specialisation of read_from_ptree<> for T
		template <typename T>
		T read_from_json_file(const boost::filesystem::path &arg_json_file ///< The JSON file to read
		                      ) {
			std::ifstream json_file_ifstream;
			open_ifstream( json_file_ifstream, arg_json_file );
			const std::string json_str{
				std::istreambuf_iterator<char>( json_file_ifstream ),
				std::istreambuf_iterator<char>()
			};
			json_file_ifstream.close();
			return from_json_string<T>( json_str );
		}

	} // namespace common
} // namespace cath

#endif
