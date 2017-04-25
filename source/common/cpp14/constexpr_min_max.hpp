/// \file
/// \brief The constexpr min()/max() class header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_CPP14_CONSTEXPR_MIN_MAX_H
#define _CATH_TOOLS_SOURCE_COMMON_CPP14_CONSTEXPR_MIN_MAX_H

#include <cstddef>

namespace cath {
	namespace common {

	/// \brief Simple constexpr min() function as workaround for GCC's min() not yet being constexpr
	///
	/// Defaults to lhs if equivalent. Only uses operator<().
	///
	/// \todo Come GCC 5.0 (or whenever libstdc++'s std::min() is constexpr) and come all Clang >= 3.5,
	///       remove this and replace all calls with calls to std::min()
	template <typename T>
	inline constexpr const T & constexpr_min(const T &arg_lhs, ///< The left-hand-side  argument
	                                         const T &arg_rhs  ///< The right-hand-side argument
	                                         ) {
		return ( arg_rhs < arg_lhs ) ? arg_rhs : arg_lhs;
	}

	/// \brief Simple constexpr max() function as workaround for GCC's max() not yet being constexpr
	///
	/// Defaults to lhs if equivalent. Only uses operator<().
	///
	/// \todo Come GCC 5.0 (or whenever libstdc++'s std::min() is constexpr) and come all Clang >= 3.5,
	///       remove this and replace all calls with calls to std::min()
	template <typename T>
	inline constexpr const T & constexpr_max(const T &arg_lhs, ///< The left-hand-side  argument
	                                         const T &arg_rhs  ///< The right-hand-side argument
	                                         ) {
		return ( arg_lhs < arg_rhs ) ? arg_rhs : arg_lhs;
	}

} // namespace common


} // namespace cath

#endif
