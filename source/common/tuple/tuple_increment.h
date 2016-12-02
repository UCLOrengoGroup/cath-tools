/// \file
/// \brief The tuple_increment header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_TUPLE_TUPLE_INCREMENT_H
#define _CATH_TOOLS_SOURCE_COMMON_TUPLE_TUPLE_INCREMENT_H

#include "common/detail/make_static_const.h"
#include "common/detail/tuple_index_sequence.h"
#include "common/tuple/is_tuple.h"

#include <tuple>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Return an incremented copy of the argument
			template <typename T>
			inline constexpr auto increment_copy(T arg ///< The value from which a copy should be taken, incremented and returned
			                                     ) {
				++arg;
				return arg;
			}

			/// \brief Implementation for tuple_increment
			template <typename Tpl, size_t... Index>
			constexpr auto tuple_increment_impl(const Tpl &arg_tuple,         ///< The tuple to be incremented
			                                    std::index_sequence<Index...> ///< An index_sequence matching the indices of Tpl
			                                    ) {
				return std::make_tuple( increment_copy( std::get<Index>( arg_tuple ) )... );
			}

			/// \brief Function object to return the result of incrementing all the members of the specified tuple
			struct tuple_increment_fn final {

				/// \brief Return the result of incrementing all the members of the specified tuple
				///
				/// \todo Tidy up the enable_if / decltype()
				template <typename Tpl, typename = std::enable_if< is_tuple< Tpl >::value > >
				constexpr auto operator()(const Tpl &arg_tuple ///< The tuple from which a copy should be taken, its members incremented and returned
				                          ) const -> decltype(
				                                     	detail::tuple_increment_impl(
				                                     		arg_tuple,
				                                     		tuple_index_sequence<Tpl>{}
				                                     	)
				                                     ) {
					return detail::tuple_increment_impl(
						arg_tuple,
						tuple_index_sequence<Tpl>{}
					);
				}

				tuple_increment_fn()                           = delete;
				tuple_increment_fn(const tuple_increment_fn &) = delete;
				void operator=(const tuple_increment_fn &)     = delete;
			};

		} // namespace detail

		MAKE_STATIC_CONST( detail::tuple_increment_fn, tuple_increment )

	} // namespace common
} // namespace cath

#endif
