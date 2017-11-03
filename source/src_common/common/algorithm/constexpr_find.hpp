/// \file
/// \brief The constexpr_find header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_CONSTEXPR_FIND_H
#define _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_CONSTEXPR_FIND_H

#include "common/cpp14/constexpr_min_max.hpp"

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

namespace cath {
	namespace common {

		namespace detail {

			/// \brief Return whether a match has been found
			///        (ie whether arg_value equals the J-th entry of the I-th entry of arg_array)
			template <size_t I, size_t J, typename T, size_t N, typename V>
			constexpr bool found_constexpr_match(const std::array<T, N> &arg_array, ///< The array to be searched
			                                     const V                &arg_value  ///< The value to be searched for
			                                     ) {
				return ( std::get<J>( std::get<I>( arg_array ) ) == arg_value );
			}

			/// \brief The recursive implementation for constexpr_find
			template <size_t I, size_t J, typename T, size_t N, typename V>
			constexpr const T & constexpr_find_impl(const std::array<T, N> &arg_array, ///< The array to be searched
			                                        const V                &arg_value  ///< The value to be searched for
			                                        ) {
				// Ensure the array isn't empty
				static_assert( N > 0, "Cannot constexpr_find() in empty std::array" );

				// * If a match is found at this position, return it
				// * Else if there are more indices left in the array, recursively try the next one
				//   (being careful to avoid even *instantiating* constexpr_find_impl() beyond the end of the array)
				// * Else fail (using a throw statement that generates compiler errors iff it's invoked)
				return found_constexpr_match<I, J>( arg_array, arg_value ) ? std::get<I>( arg_array ) :
				       ( I + 1 < N                                       ) ? constexpr_find_impl<constexpr_min(N-1, I+1), J>( arg_array, arg_value ) :
				                                                            ( throw std::logic_error( "Unable to constexpr_find() element in std::array" ), std::get<0>( arg_array ) );
			}
		} // namespace detail

		/// \brief Recursively search through an array for an entry that has arg_value in position J (default 0)
		///        (where the array's entries are of some type like pair or tuple that can be queried with std::get<J>() )
		///
		/// The important thing about this is that it is constexpr so this makes it possible to use array<pair<>> / array<tuple<>>
		/// as a compile-time map of constexpr values.
		///
		/// Tips:
		///  * Wherever there is a fixed lookup that requires values to be added for code to work, implement it as
		///    an array<pair<>> or array<tuple<>> and attempt to have the dependent code invokes a compile-time call
		///    to constexpr_find() so that the code will helpfully fail to compile on any attempt to using a missing value.
		///  * If lookups with non-constexpr values are required (eg involving string), aim to still ensure that there's at least
		///    one array of constexpr keys and a compile-time use of constexpr_find. Then add a run-time test that iterates
		///    over that array and checks all other non-constexpr lookups have values too.
		///    (Thus values missing from the constexpr lookup(s) are guaranteed to fail at compile-time and values missing from
		///     the non-constexpr lookup(s) are guaranteed to be included in the test-suite)
		///
		/// This is implemented via the recursive constexpr_find_impl()
		///
		/// \todo Come C++14 (GCC 5.0), the relaxed constexpr rules probably allow substantial
		///       simplification of the implementation code
		template <size_t J = 0, typename T, size_t N, typename V>
		constexpr const T & constexpr_find(const std::array<T, N> &arg_array, ///< The array to be searched
		                                   const V                &arg_value  ///< The value to be searched for
		                                   ) {
			return detail::constexpr_find_impl<0, J>( arg_array, arg_value );
		}

	} // namespace common
} // namespace cath

#endif
