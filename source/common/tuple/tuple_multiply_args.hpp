/// \file
/// \brief The tuple_multiply_args header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_TUPLE_TUPLE_MULTIPLY_ARGS_H
#define _CATH_TOOLS_SOURCE_COMMON_TUPLE_TUPLE_MULTIPLY_ARGS_H

#include "common/cpp17/apply.hpp"
#include "common/detail/make_static_const.hpp"
#include "common/function/multiply_args.hpp"
#include "common/type_traits/is_tuple.hpp"

#include <tuple>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Return an multiply_argsed copy of the argument
			template <typename T>
			inline constexpr auto multiply_args_copy(T arg ///< The value from which a copy should be taken, multiply_argsed and returned
			                                         ) {
				++arg;
				return arg;
			}

			/// \brief Implementation for tuple_multiply_args
			template <typename Tpl, size_t... Index>
			constexpr auto tuple_multiply_args_impl(const Tpl &arg_tuple,         ///< The tuple to be multiply_argsed
			                                        std::index_sequence<Index...> ///< An index_sequence matching the indices of Tpl
			                                        ) {
				return multiply_args( std::get<Index>( arg_tuple )... );
			}

			/// \brief Function object to return the result of multiply_argsing all the members of the specified tuple
			struct tuple_multiply_args_fn final {

				/// \brief Return the result of multiply_argsing all the members of the specified tuple
				///
				/// \todo Tidy up the enable_if / decltype()
				template <typename Tpl, typename = std::enable_if< is_tuple< Tpl >::value > >
				constexpr auto operator()(const Tpl &arg_tuple ///< The tuple from which a copy should be taken, its members multiply_argsed and returned
				                          ) const -> decltype( apply( multiply_args, arg_tuple ) ) {
					return apply( multiply_args, arg_tuple );
				}

				tuple_multiply_args_fn()                           = delete;
				tuple_multiply_args_fn(const tuple_multiply_args_fn &) = delete;
				void operator=(const tuple_multiply_args_fn &)     = delete;
			};

		} // namespace detail

		MAKE_STATIC_CONST( detail::tuple_multiply_args_fn, tuple_multiply_args )

	} // namespace common
} // namespace cath

#endif
