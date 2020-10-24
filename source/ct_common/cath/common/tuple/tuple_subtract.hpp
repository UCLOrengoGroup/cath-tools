/// \file
/// \brief The tuple_subtract header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_SUBTRACT_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_SUBTRACT_HPP

#include "cath/common/detail/make_static_const.hpp"

#include <cstddef>
#include <tuple>
#include <type_traits>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Implementation of tuple_subtract
			template <typename TplA, typename TplB, size_t... Index>
			constexpr auto tuple_subtract_impl(const TplA &prm_tuple_a,      ///< The tuple from which the other should be subtracted
			                                   const TplB &prm_tuple_b,      ///< The tuple to subtract from the other
			                                   std::index_sequence<Index...> ///< An index_sequence matching the indices of Tpl
			                                   ) {
				return std::make_tuple(
					static_cast<std::common_type_t<std::tuple_element_t<Index, TplA>, std::tuple_element_t<Index, TplA>>>(
						std::get<Index>( prm_tuple_a ) - std::get<Index>( prm_tuple_b )
					)...
				);
			}

			/// \brief Function object to return the result of element-wise subtracting one tuple from another
			struct tuple_subtract_fn final {

				/// \brief Return the result of element-wise subtracting one tuple from another
				template <typename TplA,
				          typename TplB>
				constexpr auto operator()(const TplA &prm_tuple_a, ///< The tuple from which the other should be subtracted
				                          const TplB &prm_tuple_b  ///< The tuple to subtract from the other
				                          ) const {
					static_assert(
						std::tuple_size< std::decay_t< TplA > >::value == std::tuple_size< std::decay_t< TplB > >::value,
						"tuple_subtract() can only be used on tuples of equal size"
					);
					return tuple_subtract_impl(
						prm_tuple_a,
						prm_tuple_b,
						std::make_index_sequence<std::tuple_size< std::decay_t< TplA > >::value>{}
					);
				}
			};

		} // namespace detail

		MAKE_STATIC_CONST( detail::tuple_subtract_fn, tuple_subtract )

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_SUBTRACT_HPP
