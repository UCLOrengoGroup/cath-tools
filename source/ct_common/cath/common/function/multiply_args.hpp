/// \file
/// \brief The multiply_args header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FUNCTION_MULTIPLY_ARGS_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FUNCTION_MULTIPLY_ARGS_HPP

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Implementation for multiply_args when there is one argument
			template <typename T>
			constexpr T multiply_args_impl(const T &arg ///< The argument to multiply
			                               ) {
				return arg;
			}

			/// \brief Implementation for multiply_args when there are two or more arguments
			template <typename T, typename U, typename... Vs>
			constexpr auto multiply_args_impl(const T  &   prm_1, ///< The first  argument to multiply
			                                  const U  &   prm_2, ///< The second argument to multiply
			                                  const Vs &...args   ///< All remaining arguments to multiply
			                                  ) {
				return prm_1 * multiply_args_impl<U, Vs...>( prm_2, args... );
			}

			/// \brief Function object to return the result of multiplying the arguments
			///
			/// \todo Come C++17, remove this and replace all calls with a suitable fold expression
			struct multiply_args_fn final {

				/// \brief Function operator to return the result of multiplying the arguments
				template <typename... Ts>
				constexpr auto operator()(const Ts &...args ///< The arguments to multiply
				                          ) const {
					static_assert( sizeof...(Ts) >= 1, "multiply_args must be passed at least one argument" );
					return multiply_args_impl( args... );
				}

				multiply_args_fn()                         = delete;
				multiply_args_fn(const multiply_args_fn &) = delete;
				void operator=(const multiply_args_fn &)   = delete;
			};

		} // namespace detail

		[[maybe_unused]] constexpr detail::multiply_args_fn multiply_args{};

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FUNCTION_MULTIPLY_ARGS_HPP
