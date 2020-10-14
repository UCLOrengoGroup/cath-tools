/// \file
/// \brief The tuple conjunction header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CPP17_CONJUNCTION_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CPP17_CONJUNCTION_HPP

#include <type_traits>

namespace cath {
	namespace common {

		/// \brief C++17 logical conjunction of the type traits Bs
		template <typename...> struct conjunction : std::true_type {};

		/// \brief C++17 logical conjunction of the type traits Bs
		template <typename B1> struct conjunction<B1> : B1 {};

		/// \brief C++17 logical conjunction of the type traits Bs
		template <typename B1, typename... Bn>
		struct conjunction<B1, Bn...>  : std::conditional_t< bool( B1::value ), conjunction<Bn...>, B1> {};

		// /// \brief C++17 variable template for logical conjunction of the type traits Bs
		// template <typename... B>
		// constexpr bool conjunction_v = conjunction<B...>::value;

	} // namespace common
} // namespace cath

#endif
