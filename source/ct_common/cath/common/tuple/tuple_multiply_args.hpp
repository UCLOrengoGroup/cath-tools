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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_MULTIPLY_ARGS_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_MULTIPLY_ARGS_HPP

#include "cath/common/cpp17/apply.hpp"
#include "cath/common/function/multiply_args.hpp"
#include "cath/common/type_traits/is_tuple.hpp"

#include <tuple>

namespace cath::common {
	namespace detail {

		/// \brief Function object to return the result of multiply_argsing all the members of the specified tuple
		struct tuple_multiply_args_fn final {

			/// \brief Return the result of multiply_argsing all the members of the specified tuple
			///
			/// \todo Tidy up the enable_if / decltype()
			template <typename Tpl, typename = std::enable_if< is_tuple_v< Tpl > > >
			constexpr auto operator()(const Tpl &prm_tuple ///< The tuple from which a copy should be taken, its members multiply_argsed and returned
			                          ) const -> decltype( ::cath::common::apply( multiply_args, prm_tuple ) ) {
				/// \TODO Come C++17, use ::std::apply
				return ::cath::common::apply( multiply_args, prm_tuple );
			}

			tuple_multiply_args_fn()                           = delete;
			tuple_multiply_args_fn(const tuple_multiply_args_fn &) = delete;
			void operator=(const tuple_multiply_args_fn &)     = delete;
		};

	} // namespace detail

	inline constexpr detail::tuple_multiply_args_fn tuple_multiply_args{};

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_MULTIPLY_ARGS_HPP
