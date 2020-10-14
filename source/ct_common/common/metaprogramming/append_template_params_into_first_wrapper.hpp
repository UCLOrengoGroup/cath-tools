/// \file
/// \brief The append_template_params_into_first_wrapper header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_METAPROGRAMMING_APPEND_TEMPLATE_PARAMS_INTO_FIRST_WRAPPER_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_METAPROGRAMMING_APPEND_TEMPLATE_PARAMS_INTO_FIRST_WRAPPER_HPP

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Primary template for implementation of append_template_params_into_first_wrapper_t
			template <typename... Ts>
			struct append_template_params_into_first_wrapper final {};

			/// \brief Specialisation of template for implementation of append_template_params_into_first_wrapper_t
			///        for single arguments
			template <typename T>
			struct append_template_params_into_first_wrapper<T> final {
				using type = T;
			};

			/// \brief Specialisation of template for implementation of append_template_params_into_first_wrapper_t
			///        for two or more arguments (where both must be templates)
			template <template <typename...> class Wrppr1,
			          typename... Prms1,
			          template <typename...> class Wrppr2,
			          typename... Prms2,
			          typename... Ts>
			struct append_template_params_into_first_wrapper<Wrppr1<Prms1...>, Wrppr2<Prms2...>, Ts...> final {
				using type = typename append_template_params_into_first_wrapper< Wrppr1< Prms1..., Prms2... >, Ts... >::type;
			};

		} // namespace detail

	/// \brief The type created by wrapping the specified type's template parameters in the template wrapper
	///
	/// Eg Given `tuple<int>` and `vector` this returns `vector<int>`
	template <typename... Ts>
	using append_template_params_into_first_wrapper_t = typename detail::append_template_params_into_first_wrapper<Ts...>::type;

	} // namespace common
} // namespace cath

#endif
