/// \file
/// \brief The std::array<char, N> type_aliases header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_CHAR_ARR_TYPE_ALIASES_H
#define _CATH_TOOLS_SOURCE_COMMON_CHAR_ARR_TYPE_ALIASES_H

#include <array>

namespace cath {

	/// \brief A type alias for a std::array of 2 chars
	using char_2_arr = std::array<char, 2>;

	/// \brief A type alias for a std::array of 3 chars
	using char_3_arr = std::array<char, 3>;

	/// \brief A type alias for a std::array of 4 chars
	using char_4_arr = std::array<char, 4>;

} // namespace cath

#endif
