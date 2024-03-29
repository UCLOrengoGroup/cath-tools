/// \file
/// \brief The string_view_of_char_arr header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_STRING_VIEW_OF_CHAR_ARR_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_STRING_VIEW_OF_CHAR_ARR_HPP

#include <array>
#include <string>
#include <string_view>

#include <boost/range/algorithm/find_if.hpp>

namespace cath::common {

	/// \brief Make a string_view from the specified array of chars
	template <size_t N>
	constexpr ::std::string_view string_view_of_char_arr(const std::array<char, N> &prm_char_arr ///< The array of chars from which to make the string_view
	                                                     ) {
		return { ::std::cbegin( prm_char_arr ), N };
	}

	/// \brief Make a string from the specified array of chars that
	///        may use a null to specify the end of the string
	template <size_t N>
	constexpr ::std::string_view string_view_of_null_term_char_arr(const std::array<char, N> &prm_char_arr ///< The array of chars from which to make the string
	                                                               ) {
		return {
			::std::cbegin( prm_char_arr ),
			static_cast<size_t>( std::distance(
				::std::cbegin( prm_char_arr ),
				boost::range::find_if(
					prm_char_arr,
					[] (const char &x) { return ( x == 0 ); }
				)
			) )
		};
	}

	/// \brief Make a string from the specified array of chars
	template <size_t N>
	inline std::string string_of_char_arr(const std::array<char, N> &prm_char_arr ///< The array of chars from which to make the string
	                                      ) {
		return {
			::std::cbegin( prm_char_arr ),
			::std::cend  ( prm_char_arr )
		};
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_STRING_VIEW_OF_CHAR_ARR_HPP
