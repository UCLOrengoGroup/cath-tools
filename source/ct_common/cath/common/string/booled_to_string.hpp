/// \file
/// \brief The booled_to_string header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_BOOLED_TO_STRING_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_BOOLED_TO_STRING_HPP

#include <sstream>
#include <string>

namespace cath {
	namespace common {

		/// \brief Return a string of the specified value, particularly "true"/"false" for a bool
		///
		/// This is the default template implementation that just calls std::to_string().
		template <typename T>
		inline std::string booled_to_string(const T &prm_value ///< The value to be converted to a string
		                                    ) {
			return std::to_string( prm_value );
		}

		/// \brief Specialisation to return a "true"/false" string for a bool
		///
		/// This uses std::boolalpha on a std::ostringstream.
		template <>
		inline std::string booled_to_string<bool>(const bool &prm_value ///< The bool value to be converted to a string
		                                          ) {
			std::ostringstream out_ss;
			out_ss << std::boolalpha << prm_value;
			return out_ss.str();
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_BOOLED_TO_STRING_HPP
