/// \file
/// \brief The is_template_of_type header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_IS_TEMPLATE_OF_TYPE_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_IS_TEMPLATE_OF_TYPE_HPP

#include <tuple>
#include <type_traits>

namespace cath::common {
	namespace detail {

		/// \brief Type trait for whether the first type is a template of the second template type
		///
		/// This is the generic template, which matches all types that aren't template types
		template <typename T, template <typename...> class U, typename = ::std::void_t<>>
		struct is_template_of_type : ::std::false_type {};

		/// \brief Type trait for whether the first type is a template of the second template type
		///
		/// This is the partial specialisation that matches all template types with type parameters that can be used
		/// to instantiate the specified template type. This just then checks whether the two template types
		/// are the same when instantiated with those same type parameters.
		///
		/// The use of void_t is necessary to take this specialisation out when the second (template) type
		/// can't be instantiated with the first type's template parameters. This was motivated by the case of there
		/// being more parameters than the second (template) type can accept
		/// (eg `is_template_of_type<set<int, less<int>, allocator<int>>, pair`),
		/// which otherwise causes compiler errors
		template <template <typename...> class TOut, typename... TIns, template <typename...> class U>
		struct is_template_of_type<TOut<TIns...>, U, ::std::void_t<U<TIns...>>> : ::std::is_same<U<TIns...>, TOut<TIns...>> {};

	} // namespace detail

	template <typename T, template <typename...> class U>
	constexpr bool is_template_of_type_v = detail::is_template_of_type<T, U>::value;

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TYPE_TRAITS_IS_TEMPLATE_OF_TYPE_HPP
