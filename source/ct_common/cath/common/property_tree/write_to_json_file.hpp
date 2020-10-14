/// \file
/// \brief The write_to_json_file header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_PROPERTY_TREE_WRITE_TO_JSON_FILE_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_PROPERTY_TREE_WRITE_TO_JSON_FILE_HPP

#include <boost/filesystem/path.hpp>

#include "cath/common/file/open_fstream.hpp"
#include "cath/common/json_style.hpp"
#include "cath/common/property_tree/to_json_string.hpp"

#include <fstream>

namespace cath {
	namespace common {

		/// \brief Create a JSON string of the specified value (via a ptree)
		///
		/// \tparam T must have an associated `save_to_ptree(ptree &, const T &)` non-member function
		template <typename T>
		void write_to_json_file(const boost::filesystem::path &prm_json_out_file,                  ///< TODOCUMENT
		                        const T                       &prm_val,                            ///< The value to represent in the JSON string
		                        const json_style              &prm_json_style = json_style::PRETTY ///< The style in which the JSON should be written
		                        ) {
			std::ofstream json_file_ostream;
			open_ofstream( json_file_ostream, prm_json_out_file );
			json_file_ostream << to_json_string( prm_val, prm_json_style );
			json_file_ostream << std::flush;
			json_file_ostream.close();
		}


	} // namespace common
} // namespace cath

#endif
