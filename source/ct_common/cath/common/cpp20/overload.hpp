/// \file
/// \brief The overload header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_OVERLOAD_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_OVERLOAD_HPP

#include <utility>

namespace cath::common {

	template <class... Ts>
	struct overload : Ts... {

		// TODO: Come a fix for https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80438, remove this
		constexpr explicit overload( Ts &&...params ) : Ts{ ::std::forward<Ts>( params ) }... {
		}

		using Ts::operator()...;
	};

	/// TODO: Come C++20, remove this
	template <class... Ts>
	overload( Ts... ) -> overload<Ts...>;

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_OVERLOAD_HPP
