/// \file
/// \brief The char_arr_to_string header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_STRING_CHAR_ARR_TO_STRING_H
#define _CATH_TOOLS_SOURCE_COMMON_STRING_CHAR_ARR_TO_STRING_H

#include "common/cpp14/cbegin_cend.hpp"

#include <array>
#include <string>

namespace cath {
	namespace common {

		/// \brief Return a string of the characters in the specified char array
		template <size_t N>
		inline std::string char_arr_to_string(const std::array<char, N> &arg_char_arr ///< The char array from which to build the string
		                                      ) {
			return { 
				common::cbegin( arg_char_arr ),
				common::cend  ( arg_char_arr )
			};
		}

	} // namespace common
} // namespace cath

#endif
