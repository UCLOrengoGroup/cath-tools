/// \file
/// \brief The tuple apply header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CPP17_APPLY_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CPP17_APPLY_HPP

#include <boost/core/ignore_unused.hpp> // for ignore_unused

#include "common/detail/tuple_index_sequence.hpp"

#include <cstddef>     // for size_t
#include <tuple>       // for tuple_element_t
#include <type_traits> // for decay_t
#include <utility>     // for forward, index_sequence, etc

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Implementation function for apply
			template <typename Fn, typename Tuple, std::size_t... Index>
			constexpr void apply_stepwise_impl(Fn    &&prm_fn,               ///< The function to apply to each of the tuple's elements
			                                   Tuple &&prm_tuple,            ///< The tuple of values to which the function should be applied
			                                   std::index_sequence<Index...> ///< An index_sequence matching the indices of Tpl
			                                   ) {
				auto dummy_list = {
					apply_fn_and_return_dummy_value(
						std::get<Index> ( std::forward<Tuple>( prm_tuple ) ),
						std::forward<Fn>( prm_fn )
					)...
				};
				boost::ignore_unused( dummy_list );
			}

			/// \brief Implementation function for apply
			///
			/// \todo Come C++17, remove this and switch all uses to std::apply
			template <typename Fn, typename Tpl, size_t... Index>
			constexpr decltype(auto) apply_impl(Fn  &&prm_fn,                 ///< The function to apply to the tuple's elements
			                                    Tpl &&prm_tuple,              ///< The tuple of values to which the function should be applied
			                                    std::index_sequence<Index...> ///< An index_sequence matching the indices of Tpl
			                                    ) {
				return std::forward< Fn >( prm_fn )(
					std::get< Index >( std::forward< Tpl >( prm_tuple ) )...
				);
			}

		} // namespace detail

		/// \brief Apply the function to each of the members of the tuple in order
		template <typename Fn, typename Tpl>
		constexpr void apply_stepwise(Fn  &&prm_fn,   ///< The function to apply to each of the elements
		                              Tpl &&prm_tuple ///< The tuple of values to which the function should be applied
		                              ) {
			detail::apply_stepwise_impl(
				std::forward< Fn  >( prm_fn    ),
				std::forward< Tpl >( prm_tuple ),
				detail::tuple_index_sequence<Tpl>{}
			);
		}


		/// \brief Apply the function to each of the members of the tuple in order
		///
		/// \todo Come C++17, remove this and switch all uses to std::apply
		template <typename Fn, typename Tpl>
		constexpr decltype(auto) apply(Fn  &&prm_fn,   ///< The function to apply to the tuple's elements
		                               Tpl &&prm_tuple ///< The tuple of values to which the function should be applied
		                               )  {
			return detail::apply_impl(
				std::forward< Fn  >( prm_fn    ),
				std::forward< Tpl >( prm_tuple ),
				detail::tuple_index_sequence<Tpl>{}
			);
		}

	} // namespace common
} // namespace cath

#endif
