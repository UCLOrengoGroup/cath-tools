/// \file
/// \brief The string_ref_of_char_arr header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_STRING_REF_OF_CHAR_ARR_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_STRING_REF_OF_CHAR_ARR_HPP

#include <boost/range/algorithm/find_if.hpp>
#include <boost/utility/string_ref.hpp>

#include "cath/common/cpp14/cbegin_cend.hpp"

#include <array>
#include <string>

namespace cath {
	namespace common {

		/// \brief Make a string_ref from the specified array of chars
		template <size_t N>
		inline boost::string_ref string_ref_of_char_arr(const std::array<char, N> &prm_char_arr ///< The array of chars from which to make the string_ref
		                                                ) {
			return { common::cbegin( prm_char_arr ), N };
		}

		/// \brief Make a string from the specified array of chars that
		///        may use a null to specify the end of the string
		template <size_t N>
		inline boost::string_ref string_ref_of_null_term_char_arr(const std::array<char, N> &prm_char_arr ///< The array of chars from which to make the string
		                                                          ) {
			return {
				common::cbegin( prm_char_arr ),
				static_cast<size_t>( std::distance(
					common::cbegin( prm_char_arr ),
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
				common::cbegin( prm_char_arr ),
				common::cend  ( prm_char_arr )
			};
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_STRING_REF_OF_CHAR_ARR_HPP
