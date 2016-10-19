/// \file
/// \brief The constexpr_floor header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_CONSTEXPR_FLOOR_H
#define _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_CONSTEXPR_FLOOR_H

namespace cath {
	namespace common {

		/// \brief Provide a constexpr implementation of floor whilst there's none in the standard
		///
		/// C++14 doesn't appear to offer a constexpr floor function. libstdc++'s (GCC) floor() is constexpr
		/// but libc++'s (Clang) isn't and it looks as though libc++ is better representing the standard.
		///
		/// Even if future standards don't make the standard math functions constexpr, it seems very likely
		/// that constexpr complementes will be provided. In the meantime, this provides a rough and ready
		/// constexpr floor.
		///
		/// \todo If/when the standard makes floor() constexpr or adds a new constexpr floor function,
		///       replace this with that.
		template <typename T>
		constexpr T constexpr_floor(const T &arg_value ///< The value for which the floor should be calculated
		                            ) {
			//      If ( this is a non-negative number ) then ( return the result of casting to int and back )
			// Else if ( is a whole, negative number   ) then ( return it as is )
			// Else                                           ( return the result of: negating, adding 1.0, constexpr_floor()ing and negating again)
			return ( arg_value >= 0.0 )                                              ? static_cast<T>( static_cast<int>(  arg_value ) ) :
			       ( arg_value == static_cast<T>( static_cast<int>(  arg_value ) ) ) ? arg_value                                        :
			                                                                           -constexpr_floor( static_cast<T>( 1.0 ) - arg_value );
		}

	} // namespace common
} // namespace cath

#endif
