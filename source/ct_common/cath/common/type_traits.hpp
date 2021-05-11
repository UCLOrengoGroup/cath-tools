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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_HPP

#include <type_traits>

namespace cath {
	namespace common {

		/// Type trait for the specified type with any cvref qualifiers removed
		///
		/// This is like ::std::decay but it doesn't also decay arrays to pointers etc
		///
		/// TODO: Come C++20, retire this in favour of ::std::remove_cvref
		template <typename T>
		struct remove_cvref {
			typedef ::std::remove_cv_t<::std::remove_reference_t<T>> type;
		};

		/// Type alias for the specified type with any cvref qualifiers removed
		///
		/// This is like ::std::decay but it doesn't also decay arrays to pointers etc
		///
		/// TODO: Come C++20, retire this in favour of ::std::remove_cvref_t
		template <typename T>
		using remove_cvref_t = typename remove_cvref<T>::type;

		/// Whether the two specified types are the same, modulo cv-ref qualifiers
		///
		/// TODO: Come C++17, add variable template is_same_modulo_cvref_v as:
		///           template <typename T, typename U>
		///           inline constexpr bool is_same_modulo_cvref_v = is_same_modulo_cvref<T, U>::value;
		///       and replace most uses of this with uses of that.
		template <typename T, typename U>
		struct is_same_modulo_cvref {
			static constexpr bool value = ::std::is_same<remove_cvref_t<T>, remove_cvref_t<U>>::value;
		};

		template <typename T, typename U>
		inline constexpr bool is_same_modulo_cvref_v = is_same_modulo_cvref<T, U>::value;

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_IS_TUPLE_HPP
