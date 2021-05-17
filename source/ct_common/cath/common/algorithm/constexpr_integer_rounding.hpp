/// \file
/// \brief The constexpr_modulo_fns header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_CONSTEXPR_INTEGER_ROUNDING_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_CONSTEXPR_INTEGER_ROUNDING_HPP

#include <type_traits>
	
namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		template <typename T>
		inline constexpr T round_div_up(T a, ///< TODOCUMENT
		                                T b  ///< TODOCUMENT
		                                ) {
			static_assert( std::is_integral_v<T>, "round_div_down() must be performed on an integral type" );
			return ( a < 0 || b < 0 ) ? throw("round_div_down() currently only implemented for non-negative values")
			                          : ( a / b ) + ( ( a % b != 0 ) ? 1 : 0);
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline constexpr T round_down_mod(T a, ///< TODOCUMENT
		                                  T b  ///< TODOCUMENT
		                                  ) {
			static_assert( std::is_integral_v<T>, "round_down_mod() must be performed on an integral type" );
			return ( a < 0 || b < 0 ) ? throw("round_down_mod() currently only implemented for non-negative values")
			                          : b * ( a / b );
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline constexpr T round_up_mod(T a, ///< TODOCUMENT
		                                T b  ///< TODOCUMENT
		                                ) {
			static_assert( std::is_integral_v<T>, "round_up_mod() must be performed on an integral type" );
			return ( a < 0 || b < 0 ) ? throw("round_up_mod() currently only implemented for non-negative values")
			                          : b * round_div_up( a, b );
		}
	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_CONSTEXPR_INTEGER_ROUNDING_HPP
