/// \file
/// \brief The constexpr_is_uniq class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_CONSTEXPR_IS_UNIQ_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_CONSTEXPR_IS_UNIQ_HPP

#include <array>
#include <cstddef>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief TODOCUMENT
			template <size_t A, size_t N>
			constexpr size_t constexpr_is_uniq_impl_next_index() {
				return ( A + 1 < N ) ? A + 1 : A;
			}

			/// \brief TODOCUMENT
			template <size_t I = 0,size_t J = 0, typename T, size_t N>
			constexpr bool constexpr_is_uniq_impl_impl(const std::array<T, N> &prm_array ///< The array to be searched
			                                           ) {
				return (
					( I == J || prm_array[ I ] != prm_array [ J ] )
					&&
					(
						( J + 1 >= N )
						||
						constexpr_is_uniq_impl_impl< I, constexpr_is_uniq_impl_next_index<J, N>() >( prm_array )
					)
				);
			}

			/// \brief TODOCUMENT
			template <size_t I = 0, typename T, size_t N>
			constexpr bool constexpr_is_uniq_impl(const std::array<T, N> &prm_array ///< The array to be searched
			                                      ) {
				return (
					constexpr_is_uniq_impl_impl<I, 0>( prm_array )
					&&
					(
						( I + 1 >= N )
						||
						constexpr_is_uniq_impl<constexpr_is_uniq_impl_next_index<I, N>()>( prm_array )
					)
				);
			}

		} // namespace detail

		/// \brief TODOCUMENT
		template <typename T, size_t N>
		constexpr bool constexpr_is_uniq(const std::array<T, N> &prm_array ///< The array to be searched
		                                 ) {
			return detail::constexpr_is_uniq_impl( prm_array );
		}


//
//			/// \brief TODOCUMENT
//			template <size_t I, size_t J, size_t N>
//			constexpr size_t constexpr_is_uniq_impl_next_i() {
//				return ( J + 1 >= N && I + 1 <  N )  ?  I + 1
//				                                     :  I;
//			}
//
//			/// \brief TODOCUMENT
//			template <size_t I, size_t J, size_t N>
//			constexpr size_t constexpr_is_uniq_impl_next_j() {
//				return ( I + 1 >= N && J + 1 >= N )  ?  J      :
//				       (               J + 1 >= N )  ?  0      :
//				                                        J + 1;
//			}
//
//			/// \brief TODOCUMENT
//			template <size_t I = 0, size_t J = 0, typename T, size_t N>
//			constexpr bool constexpr_is_uniq_impl(const std::array<T, N> &prm_array ///< The array to be searched
//			                                      ) {
//				return ( I != J && prm_array[ I ] == prm_array [ J ] ) ? false :
//				       ( I + 1 >= N && J + 1 >= N                    ) ? true  :
//				                                                         constexpr_is_uniq_impl< constexpr_is_uniq_impl_next_i<I, J, N>(),
//				                                                                                 constexpr_is_uniq_impl_next_j<I, J, N>() >( prm_array );
//			}

	} // namespace common
} // namespace cath

#endif
