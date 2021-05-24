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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_IS_TUPLE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_IS_TUPLE_HPP

#include <tuple>

#include "cath/common/type_traits.hpp"
#include "cath/common/type_traits/is_template_of_type.hpp"

namespace cath::common {

	/// \brief A type_trait for checking whether T is exactly a std::tuple<> of zero or more types
	template <typename T>
	constexpr bool is_tuple_v = is_template_of_type_v<T, ::std::tuple>;

	/// \brief A type_trait for checking whether the decayed version of T is a std::tuple<> of zero or more types
	template <typename T>
	constexpr bool is_tuple_mod_cvref_v = is_template_of_type_v<common::remove_cvref_t<T>, ::std::tuple>;

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_IS_TUPLE_HPP
