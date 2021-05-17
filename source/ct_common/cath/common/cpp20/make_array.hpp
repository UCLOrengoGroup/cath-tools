/// \file
/// \brief The tuple make_array header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_MAKE_ARRAY_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_MAKE_ARRAY_HPP

#include <array>

#include "cath/common/cpp17/conjunction.hpp"
#include "cath/common/cpp17/negation.hpp"
#include "cath/common/type_traits.hpp"

namespace cath {
	namespace common {
		namespace detail {

			template <typename> struct is_ref_wrapper : std::false_type {};

			template <typename T> struct is_ref_wrapper<std::reference_wrapper<T>> : std::true_type {};

			template <typename T>
			using not_ref_wrapper = negation<is_ref_wrapper<common::remove_cvref_t<T>>>;


			template <typename D, typename...> struct return_type_helper { using type = D; };


			template <typename... Types>
			struct return_type_helper<void, Types...> : std::common_type<Types...> {
				static_assert(::std::conjunction_v<not_ref_wrapper<Types>...>,
				"Types cannot contain reference_wrappers when D is void");
			};

			template <typename D, typename... Types>
			using return_type = std::array<typename return_type_helper<D, Types...>::type, sizeof...(Types)>;
		}

		/// \brief C++20? factory function for arrays
		template <typename D = void, typename... Types>
		constexpr detail::return_type<D, Types...> make_array(Types&&... t) {
			return { { std::forward<Types>(t)... } };
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_MAKE_ARRAY_HPP
