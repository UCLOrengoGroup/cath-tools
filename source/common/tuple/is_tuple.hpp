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

#ifndef _CATH_TOOLS_SOURCE_COMMON_TUPLE_IS_TUPLE_H
#define _CATH_TOOLS_SOURCE_COMMON_TUPLE_IS_TUPLE_H

#include <tuple>
#include <type_traits>

namespace cath {
	namespace common {

		/// \brief A type_trait for checking whether T is exactly a std::tuple<> of zero or more types
		///
		/// This is the general case that inherits from false_type
		template <typename T>       struct is_tuple                       final : std::false_type {};

		/// \brief A type_trait for checking whether T is exactly a std::tuple<> of zero or more types
		///
		/// This is the specialised case that inherits from true_type when T matches
		/// std::tuple<ARGS...> for some ...ARGS
		template <typename... Args> struct is_tuple< std::tuple<Args...> > final : std::true_type  {};

		/// \brief A type_trait for checking whether the decayed version of T is a std::tuple<> of zero or more types
		template <typename T>
		using is_tuple_after_decay = is_tuple< std::decay_t<T> >;

	} // namespace common
} // namespace cath

#endif
