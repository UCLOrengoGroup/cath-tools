/// \file
/// \brief The transform_tuple() header

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

#ifndef TRANSFORM_TUPLE_H_INCLUDED
#define TRANSFORM_TUPLE_H_INCLUDED

#include <tuple>
#include <type_traits>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief TODOCUMENT
			template <typename... Ts, typename F, size_t ...S>
			constexpr auto transform_tuple_impl(std::tuple<Ts...>         &arg_tuple,  ///< TODOCUMENT
			                                    F                          arg_fn,     ///< TODOCUMENT
			                                    std::index_sequence<S...>  /*arg_indices*/ ///< TODOCUMENT
			                                    ) {
				return arg_fn( std::get<S>( arg_tuple )... );
			}

			/// \brief TODOCUMENT
			template <typename... Ts, typename F, size_t ...S>
			constexpr auto transform_tuple_impl(const std::tuple<Ts...>   &arg_tuple,  ///< TODOCUMENT
			                                    F                          arg_fn,     ///< TODOCUMENT
			                                    std::index_sequence<S...>  /*arg_indices*/ ///< TODOCUMENT
			                                    ) {
				return arg_fn( std::get<S>( arg_tuple )... );
			}
		}

		/// \brief TODOCUMENT
		template <typename... Ts, typename F>
		constexpr auto transform_tuple(std::tuple<Ts...> &arg_tuple, ///< TODOCUMENT
		                               F                  arg_fn     ///< TODOCUMENT
		                               ) {
			return detail::transform_tuple_impl( arg_tuple, arg_fn, std::index_sequence_for<Ts...>() );
		}

		/// \brief TODOCUMENT
		template <typename... Ts, typename F>
		constexpr auto transform_tuple(const std::tuple<Ts...> &arg_tuple, ///< TODOCUMENT
		                               F                        arg_fn     ///< TODOCUMENT
		                               ) {
			return detail::transform_tuple_impl( arg_tuple, arg_fn, std::index_sequence_for<Ts...>() );
		}

//		/// \brief TODOCUMENT
//		struct lvalue_ref_tuple_maker final {
//			template <typename... RNGs>
//			std::tuple<std::add_lvalue_reference_t<RNGs>...> operator()(RNGs &... arg_ranges
//							                                            ) {
//				return { arg_ranges... };
//			}
//		};

	}
}

#endif
