/// \file
/// \brief The variadic_and header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_VARIADIC_AND_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_VARIADIC_AND_HPP

namespace cath::common {
	namespace detail {

		/// \brief One argument implementation for and-ing together all arguments
		constexpr bool variadic_and_impl(const bool &arg ///< The argument to and together
		                                 ) {
			return arg;
		}

		/// \brief Two+ argument implementation for and-ing together all arguments
		template <typename T, typename U, typename... Vs>
		constexpr bool variadic_and_impl(const T  &prm_1,  ///< The first argument to and together
		                                 const U  &prm_2,  ///< The second argument to and together
		                                 const Vs &...args ///< Any remaining arguments to and together
		                                 ) {
			return prm_1 && variadic_and_impl( prm_2, args... );
		}

		/// \brief Function object to logically "and" all the arguments
		///
		/// This isn't very smart about the arguments being bools - just use bools
		///
		/// \todo Come C++17, remove this and replace all uses with a suitable fold expression
		struct variadic_and_fn final {

			/// \brief Logically "and" all the arguments
			template <typename... Ts>
			constexpr bool operator()(const Ts &...args ///< The arguments to logically "and"
			                          ) const {
				return variadic_and_impl( args... );
			}

			variadic_and_fn()                        = delete;
			variadic_and_fn(const variadic_and_fn &) = delete;
			void operator=(const variadic_and_fn &)  = delete;
		};

	} // namespace detail

	inline constexpr detail::variadic_and_fn variadic_and{};

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_VARIADIC_AND_HPP
