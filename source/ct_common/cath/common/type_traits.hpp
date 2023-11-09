/// \file
/// \brief The is_tuple header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_HPP

#include <type_traits>

namespace cath::common {

	/// Type alias for the specified type with any cvref qualifiers removed
	///
	/// This is like ::std::decay but it doesn't also decay arrays to pointers etc
	///
	/// TODO: Come C++20, retire this in favour of ::std::remove_cvref_t
	template <typename T>
	using remove_cvref_t = ::std::remove_cv_t<::std::remove_reference_t<T>>;

	template <typename T, typename U>
	constexpr bool is_same_modulo_cvref_v = ::std::is_same_v<remove_cvref_t<T>, remove_cvref_t<U>>;

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_HPP
