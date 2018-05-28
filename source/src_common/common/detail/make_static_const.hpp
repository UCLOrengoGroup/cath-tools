/// \file
/// \brief The make_static_const header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_DETAIL_MAKE_STATIC_CONST_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_DETAIL_MAKE_STATIC_CONST_HPP

namespace cath {
	namespace detail {

		/// \brief Simple trick to avoid ODR violations
		///
		/// Hat-tip: Eric Niebler, eg see https://github.com/ericniebler/range-v3
		template <typename T>
		struct make_static_const {
			static constexpr T value{};
		};

		/// \brief Definition of the make_static_const value
		template<typename T>
		constexpr T make_static_const<T>::value;

	} // namespace detail
} // namespace cath

#define MAKE_STATIC_CONST( type, name ) inline namespace { constexpr auto const &name = cath::detail::make_static_const<type>::value; }

#endif
