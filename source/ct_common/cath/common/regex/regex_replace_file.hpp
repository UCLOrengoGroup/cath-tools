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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_REGEX_REGEX_REPLACE_FILE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_REGEX_REGEX_REPLACE_FILE_HPP

#include <boost/filesystem/path.hpp>

#include "cath/common/file/slurp.hpp"
#include "cath/common/file/spew.hpp"

#include <regex>

namespace cath {
	namespace common {

		/// \brief Perform an in-place find/replace on the specified file
		inline void regex_replace_file(const boost::filesystem::path &prm_file,  ///< The file to modify
		                               const std::regex              &prm_regex, ///< The find regex, to be passed to std::regex_replace()
		                               const std::string             &prm_string ///< The replace string, to be passed to std::regex_replace()
		                               ) {
			spew(
				prm_file,
				std::regex_replace(
					slurp( prm_file ),
					prm_regex,
					prm_string
				)
			);
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_REGEX_REGEX_REPLACE_FILE_HPP
