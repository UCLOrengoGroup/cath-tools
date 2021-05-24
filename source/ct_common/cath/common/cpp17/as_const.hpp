/// \file
/// \brief The as_const header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP17_AS_CONST_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP17_AS_CONST_HPP

#include <type_traits>

namespace cath::common {

	/// \brief C++17 helper to return a const lvalue ref to the specified lvalue argument
	template <class T>
	constexpr std::add_const_t<T> & as_const(T& t) noexcept {
		return t;
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP17_AS_CONST_HPP
