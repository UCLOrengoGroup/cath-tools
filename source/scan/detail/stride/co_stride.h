/// \file
/// \brief The co_stride header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
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

#ifndef CO_STRIDE_H_INCLUDED
#define CO_STRIDE_H_INCLUDED

#include <boost/optional.hpp>

#include "common/algorithm/constexpr_modulo_fns.h"

#include <type_traits>
#include <utility>

namespace cath {
	namespace scan {
		namespace detail {
			namespace detail {

				/// \brief TODOCUMENT
				template <typename T>
				inline constexpr T stride_neighbour_index_of_centre(const T &arg_co_stride ///< TODOCUMENT
				                                                    ) {
					static_assert( std::is_unsigned<T>::value, "stride_neighbour_index_of_centre() must be performed on an unsigned integral type" );
					return arg_co_stride / 2;
				}

				/// \brief TODOCUMENT
				///
				/// \returns TODOCUMENT
				///
				/// \todo When there is a constexpr std::optional<> available, use that as the return type instead
				template <typename T>
				inline constexpr std::pair<bool, T> entry_index_of_stride_neighbour_index_impl(const T &arg_stride_index,       ///< TODOCUMENT
				                                                                               const T &arg_co_stride,          ///< TODOCUMENT
				                                                                               const T &arg_centre_entry_index, ///< TODOCUMENT
				                                                                               const T &arg_num_entries         ///< TODOCUMENT
				                                                                               ) {
					static_assert( std::is_unsigned<T>::value, "entry_index_of_stride_neighbour_index_impl() must be performed on an unsigned integral type" );
					return               ( arg_centre_entry_index + arg_stride_index <  detail::stride_neighbour_index_of_centre( arg_co_stride )                   ) ? std::pair<bool, T>( false, 0 ) :
					                     ( arg_centre_entry_index + arg_stride_index >= detail::stride_neighbour_index_of_centre( arg_co_stride ) + arg_num_entries ) ? std::pair<bool, T>( false, 0 ) :
					 std::make_pair( true, arg_centre_entry_index + arg_stride_index -  detail::stride_neighbour_index_of_centre( arg_co_stride ) );
				}

			}

			/// \brief TODOCUMENT
			template <typename T>
			inline constexpr T co_stride(const T &arg_stride_a, ///< The stride TODOCUMENT
			                             const T &arg_stride_b  ///< The stride TODOCUMENT
			                             ) {
				static_assert( std::is_unsigned<T>::value, "co_stride() must be performed on an unsigned integral type" );
				return common::constexpr_lcm( arg_stride_a + 1, arg_stride_b + 1 ) - 1;
			}

			/// \brief TODOCUMENT
			template <typename T>
			inline constexpr T entry_index_of_stride_rep(const T &arg_entry_index, ///< TODOCUMENT
			                                             const T &arg_co_stride    ///< TODOCUMENT
			                                             ) {
				return ( arg_co_stride + 1 ) * ( ( arg_entry_index + detail::stride_neighbour_index_of_centre( arg_co_stride ) ) / ( arg_co_stride + 1 ) );
			}

			/// \brief TODOCUMENT
			template <typename T>
			inline constexpr T num_in_stride_neighbour_range(const T &arg_co_stride ///< TODOCUMENT
			                                                 ) {
				static_assert( std::is_unsigned<T>::value, "num_stride_neighbour_range() must be performed on an unsigned integral type" );
				return arg_co_stride + 1;
			}

			/// \brief TODOCUMENT
			///
			/// \returns TODOCUMENT
			template <typename T>
			inline boost::optional<T> entry_index_of_stride_neighbour_index(const T &arg_stride_index,       ///< TODOCUMENT
			                                                                const T &arg_co_stride,          ///< TODOCUMENT
			                                                                const T &arg_centre_entry_index, ///< TODOCUMENT
			                                                                const T &arg_num_entries         ///< TODOCUMENT
			                                                                ) {
				static_assert( std::is_unsigned<T>::value, "entry_index_of_stride_neighbour_index() must be performed on an unsigned integral type" );
				const auto result = detail::entry_index_of_stride_neighbour_index_impl( arg_stride_index, arg_co_stride, arg_centre_entry_index, arg_num_entries );
				return result.first ? boost::make_optional( result.second ) : boost::none;
			}



		}
	}
}

#endif
