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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CHAR_ARR_TYPE_ALIASES_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CHAR_ARR_TYPE_ALIASES_HPP

#include <boost/optional/optional_fwd.hpp>

#include <array>

namespace cath {

	namespace detail {

		/// \brief Implementation of make_char_arr() with the indices of the integers as template parameters
		template <size_t Length, size_t... Indices>
		constexpr std::array<char, Length - 1> make_char_arr_impl(const char ( &prm_cstring )[ Length ],    ///< The C-string from which the std::array should be built
		                                                          std::integer_sequence<size_t, Indices...> ///< integer_sequence encoding the indices of the individual characters
		                                                          ) {
			return { { prm_cstring[ Indices ]... } };
		}

	} // namespace detail

	/// \brief Make a std::array of the characters in the specified C-string
	template <size_t Length>
	constexpr auto make_char_arr(const char ( &prm_cstring )[ Length ] ///< The C-string from which the std::array should be built
	                             ) {
		return detail::make_char_arr_impl( prm_cstring, std::make_index_sequence<Length - 1>{} );
	}

	/// \brief A type alias for a std::array of 2 chars
	using char_2_arr = std::array<char, 2>;

	/// \brief A type alias for a std::array of 3 chars
	using char_3_arr = std::array<char, 3>;

	/// \brief A type alias for a std::array of 4 chars
	using char_4_arr = std::array<char, 4>;

	/// \brief A type alias for an optional std::array of 3 chars
	using char_3_arr_opt = boost::optional<std::array<char, 3>>;

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CHAR_ARR_TYPE_ALIASES_HPP
