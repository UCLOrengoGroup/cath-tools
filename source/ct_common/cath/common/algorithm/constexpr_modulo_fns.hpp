/// \file
/// \brief The constexpr_modulo_fns header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_CONSTEXPR_MODULO_FNS_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_CONSTEXPR_MODULO_FNS_HPP

#include "cath/common/algorithm/constexpr_integer_rounding.hpp"
#include "cath/common/cpp14/constexpr_min_max.hpp"
#include "cath/common/type_aliases.hpp"

#include <type_traits>
	
namespace cath {
	namespace common {
		template <typename T>
		inline constexpr T constexpr_gcd(T, T);
		template <typename T>
		inline constexpr T constexpr_lcm(T, T);

		namespace detail {

			/// \brief TODOCUMENT
			template <typename T, typename U>
			inline constexpr std::pair<U, T> flip_pair(const std::pair<T, U> &prm_pair ///< TODOCUMENT
			                                           ) {
				return std::make_pair( prm_pair.second, prm_pair.first );
			}

			/// \brief TODOCUMENT
			///
			/// \pre prm_r_i >= prm_r_i_plus_one
			template <typename T, typename U>
			inline constexpr std::pair<U, U> extended_euclid_algo_impl(const T &prm_r_i,          ///< TODOCUMENT
			                                                           const U &prm_s_i,          ///< TODOCUMENT
			                                                           const U &prm_t_i,          ///< TODOCUMENT
			                                                           const T &prm_r_i_plus_one, ///< TODOCUMENT
			                                                           const U &prm_s_i_plus_one, ///< TODOCUMENT
			                                                           const U &prm_t_i_plus_one  ///< TODOCUMENT
			                                                           ) {
				static_assert( std::is_unsigned<T>::value, "constexpr_gcd() must be performed on an unsigned integral type" );
				return ( prm_r_i_plus_one == 0 ) ? throw("Cannot perform extended_euclid_algo_impl with a 0") :
				       ( prm_r_i_plus_one == 1 ) ? std::make_pair( prm_s_i_plus_one, prm_t_i_plus_one ) :
				                                   extended_euclid_algo_impl(
				                                   	prm_r_i_plus_one,
				                                   	prm_s_i_plus_one,
				                                   	prm_t_i_plus_one,
				                                   	static_cast<T>( prm_r_i % prm_r_i_plus_one ),
				                                   	static_cast<U>( prm_s_i - ( prm_s_i_plus_one * static_cast<U>( prm_r_i / prm_r_i_plus_one) ) ),
				                                   	static_cast<U>( prm_t_i - ( prm_t_i_plus_one * static_cast<U>( prm_r_i / prm_r_i_plus_one) ) )
				                                   );
			}

			/// \brief TODOCUMENT
			template <typename T>
			inline constexpr auto extended_euclid_algo(const T &prm_a, ///< TODOCUMENT
			                                           const T &prm_b  ///< TODOCUMENT
			                                           ) {
				return ( prm_a >= prm_b )
				? extended_euclid_algo_impl(
					prm_a, static_cast<std::make_signed_t<T>>( 1 ), static_cast<std::make_signed_t<T>>( 0 ),
					prm_b, static_cast<std::make_signed_t<T>>( 0 ), static_cast<std::make_signed_t<T>>( 1 )
				)
				: flip_pair( extended_euclid_algo_impl(
					prm_b, static_cast<std::make_signed_t<T>>( 1 ), static_cast<std::make_signed_t<T>>( 0 ),
					prm_a, static_cast<std::make_signed_t<T>>( 0 ), static_cast<std::make_signed_t<T>>( 1 )
				) );
			}

			/// \brief TODOCUMENT
			template <typename T>
			inline constexpr auto extended_euclid_algo_products(const T &prm_a, ///< TODOCUMENT
			                                                    const T &prm_b  ///< TODOCUMENT
			                                                    ) {
				return std::make_pair(
					detail::extended_euclid_algo( prm_a, prm_b ).first  * static_cast<std::make_signed_t<T>>( prm_a ),
					detail::extended_euclid_algo( prm_a, prm_b ).second * static_cast<std::make_signed_t<T>>( prm_b )
				);
			}

			/// \brief TODOCUMENT
			template <typename T>
			inline constexpr std::make_signed_t<T> chinese_remainder_coprime_pair_num(const T &prm_index_a, ///< TODOCUMENT
			                                                                          const T &prm_index_b, ///< TODOCUMENT
			                                                                          const T &prm_mod_a,   ///< TODOCUMENT
			                                                                          const T &prm_mod_b    ///< TODOCUMENT
			                                                                          ) {
				return (
					( detail::extended_euclid_algo_products( prm_mod_a, prm_mod_b ).first  * static_cast<std::make_signed_t<T>>( prm_index_a ) )
					+
					( detail::extended_euclid_algo_products( prm_mod_a, prm_mod_b ).second * static_cast<std::make_signed_t<T>>( prm_index_b ) )
				) % static_cast<std::make_signed_t<T>>( constexpr_lcm( prm_mod_a, prm_mod_b ) );
			}

			/// \brief TODOCUMENT
			template <typename T>
			inline constexpr T chinese_remainder_coprime_pair_num_above(const T &prm_index_a, ///< TODOCUMENT
			                                                            const T &prm_index_b, ///< TODOCUMENT
			                                                            const T &prm_mod_a,   ///< TODOCUMENT
			                                                            const T &prm_mod_b    ///< TODOCUMENT
			                                                            ) {
				return
					(
						static_cast<std::make_signed_t<T>>( prm_index_a + prm_index_b )
						>
						chinese_remainder_coprime_pair_num( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b )
					)
					? static_cast<T>(
						chinese_remainder_coprime_pair_num( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b )
						+
						(
							round_up_mod(
								static_cast<std::make_signed_t<T>>( prm_index_a + prm_index_b )
								-
								chinese_remainder_coprime_pair_num( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b ),
								static_cast<std::make_signed_t<T>>( constexpr_lcm( prm_mod_a, prm_mod_b ) )
							)
						)
					)
					: static_cast<T>(
						chinese_remainder_coprime_pair_num( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b )
					);
			}

//			/// \brief TODOCUMENT
//			template <typename T>
//			inline constexpr T chinese_remainder_coprime_pair_num_above(const T &prm_index_a, ///< TODOCUMENT
//			                                                            const T &prm_index_b, ///< TODOCUMENT
//			                                                            const T &prm_mod_a,   ///< TODOCUMENT
//			                                                            const T &prm_mod_b    ///< TODOCUMENT
//			                                                            ) {
//				return
//					(
//						static_cast<std::make_signed_t<T>>( prm_index_a + prm_index_b )
//						>
//						chinese_remainder_coprime_pair_num( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b )
//					)
//					? static_cast<T>(
//						chinese_remainder_coprime_pair_num( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b )
//						+
//						(
//							static_cast<std::make_signed_t<T>>( constexpr_lcm( prm_mod_a, prm_mod_b ) )
//							*
//							(
//								static_cast<std::make_signed_t<T>>( 1 )
//								+
//								(
//									(
//										static_cast<std::make_signed_t<T>>( prm_index_a + prm_index_b )
//										-
//										chinese_remainder_coprime_pair_num( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b )
//									)
//									/
//									static_cast<std::make_signed_t<T>>( constexpr_lcm( prm_mod_a, prm_mod_b ) )
//								)
//							)
//						)
//					)
//					: static_cast<T>(
//						chinese_remainder_coprime_pair_num( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b )
//					);
//			}

			/// \brief TODOCUMENT
			template <typename T>
			inline constexpr std::pair<T, T> chinese_remainder_coprime_pair(const T &prm_index_a, ///< TODOCUMENT
			                                                                const T &prm_index_b, ///< TODOCUMENT
			                                                                const T &prm_mod_a,   ///< TODOCUMENT
			                                                                const T &prm_mod_b    ///< TODOCUMENT
			                                                                ) {
				return std::make_pair(
					chinese_remainder_coprime_pair_num_above( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b ) - prm_index_b,
					chinese_remainder_coprime_pair_num_above( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b ) - prm_index_a
				);
			}
		} // namespace detail

		/// \brief TODOCUMENT
		///
		/// \todo Come C++17, switch to std::gdd()
		template <typename T>
		inline constexpr T constexpr_gcd(T a, ///< TODOCUMENT
		                                 T b  ///< TODOCUMENT
		                                 ) {
			static_assert( std::is_unsigned<T>::value, "constexpr_gcd() must be performed on an unsigned integral type" );
			return ( b == 0 ) ? a
			                  : constexpr_gcd( b, a % b );
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come C++17, switch to std::lcm()
		template <typename T>
		inline constexpr T constexpr_lcm(T a, ///< TODOCUMENT
		                                 T b  ///< TODOCUMENT
		                                 ) {
			static_assert( std::is_unsigned<T>::value, "constexpr_lcm() must be performed on an unsigned integral type" );
			return ( a != 0 && b != 0 ) ? ( a / constexpr_gcd( a, b ) ) * b
			                            : 0;
		}

		// template <typename T> class TD;

		/// \brief TODOCUMENT
		template <typename T>
		inline constexpr std::pair<T, T> chinese_remainder_coprime_pair(const T &prm_index_a, ///< TODOCUMENT
		                                                                const T &prm_index_b, ///< TODOCUMENT
		                                                                const T &prm_mod_a,   ///< TODOCUMENT
		                                                                const T &prm_mod_b    ///< TODOCUMENT
		                                                                ) {
			return detail::chinese_remainder_coprime_pair( prm_index_a, prm_index_b, prm_mod_a, prm_mod_b );
		}

	} // namespace common
} // namespace cath

#endif
