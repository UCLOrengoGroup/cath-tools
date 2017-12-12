/// \file
/// \brief The regex_replace_file header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_REGEX_REGEX_REPLACE_FILE_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_REGEX_REGEX_REPLACE_FILE_H

#include <boost/filesystem/path.hpp>

#include "common/file/slurp.hpp"
#include "common/file/spew.hpp"

#include <regex>

namespace cath {
	namespace common {

		/// \brief Perform an in-place find/replace on the specified file
		void regex_replace_file(const boost::filesystem::path &arg_file,  ///< The file to modify
		                        const std::regex              &arg_regex, ///< The find regex, to be passed to std::regex_replace()
		                        const std::string             &arg_string ///< The replace string, to be passed to std::regex_replace()
		                        ) {
			spew(
				arg_file,
				std::regex_replace(
					slurp( arg_file ),
					arg_regex,
					arg_string
				)
			);
		}

	} // namespace common
} // namespace cath

#endif
