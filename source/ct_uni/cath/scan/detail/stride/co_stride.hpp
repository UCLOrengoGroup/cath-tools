/// \file
/// \brief The co_stride header

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_STRIDE_CO_STRIDE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_STRIDE_CO_STRIDE_HPP

#include <numeric>
#include <optional>
#include <type_traits>
#include <utility>

#include "cath/common/algorithm/constexpr_modulo_fns.hpp"

namespace cath::scan::detail {

	namespace detail {

		/// \brief TODOCUMENT
		template <typename T>
		inline constexpr T stride_neighbour_index_of_centre(const T &prm_co_stride ///< TODOCUMENT
		                                                    ) {
			static_assert( std::is_unsigned_v<T>, "stride_neighbour_index_of_centre() must be performed on an unsigned integral type" );
			return prm_co_stride / 2;
		}

		/// \brief TODOCUMENT
		///
		/// \returns TODOCUMENT
		///
		/// \todo When there is a constexpr std::optional<> available, use that as the return type instead
		template <typename T>
		inline constexpr std::pair<bool, T> entry_index_of_stride_neighbour_index_impl(const T &prm_stride_index,       ///< TODOCUMENT
		                                                                               const T &prm_co_stride,          ///< TODOCUMENT
		                                                                               const T &prm_centre_entry_index, ///< TODOCUMENT
		                                                                               const T &prm_num_entries         ///< TODOCUMENT
		                                                                               ) {
			static_assert( std::is_unsigned_v<T>, "entry_index_of_stride_neighbour_index_impl() must be performed on an unsigned integral type" );
			return               ( prm_centre_entry_index + prm_stride_index <  detail::stride_neighbour_index_of_centre( prm_co_stride )                   ) ? std::pair<bool, T>( false, 0 ) :
			                     ( prm_centre_entry_index + prm_stride_index >= detail::stride_neighbour_index_of_centre( prm_co_stride ) + prm_num_entries ) ? std::pair<bool, T>( false, 0 ) :
			 std::make_pair( true, prm_centre_entry_index + prm_stride_index -  detail::stride_neighbour_index_of_centre( prm_co_stride ) );
		}

	} // namespace detail

	/// \brief TODOCUMENT
	template <typename T>
	inline constexpr T co_stride(const T &prm_stride_a, ///< The stride TODOCUMENT
	                             const T &prm_stride_b  ///< The stride TODOCUMENT
	                             ) {
		static_assert( std::is_unsigned_v<T>, "co_stride() must be performed on an unsigned integral type" );
		return ::std::lcm( prm_stride_a + 1, prm_stride_b + 1 ) - 1;
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline constexpr T entry_index_of_stride_rep(const T &prm_entry_index, ///< TODOCUMENT
	                                             const T &prm_co_stride    ///< TODOCUMENT
	                                             ) {
		return ( prm_co_stride + 1 ) * ( ( prm_entry_index + detail::stride_neighbour_index_of_centre( prm_co_stride ) ) / ( prm_co_stride + 1 ) );
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline constexpr T num_in_stride_neighbour_range(const T &prm_co_stride ///< TODOCUMENT
	                                                 ) {
		static_assert( std::is_unsigned_v<T>, "num_stride_neighbour_range() must be performed on an unsigned integral type" );
		return prm_co_stride + 1;
	}

	/// \brief TODOCUMENT
	///
	/// \returns TODOCUMENT
	template <typename T>
	inline ::std::optional<T> entry_index_of_stride_neighbour_index(const T &prm_stride_index,       ///< TODOCUMENT
	                                                                const T &prm_co_stride,          ///< TODOCUMENT
	                                                                const T &prm_centre_entry_index, ///< TODOCUMENT
	                                                                const T &prm_num_entries         ///< TODOCUMENT
	                                                                ) {
		static_assert( std::is_unsigned_v<T>, "entry_index_of_stride_neighbour_index() must be performed on an unsigned integral type" );
		const auto result = detail::entry_index_of_stride_neighbour_index_impl( prm_stride_index, prm_co_stride, prm_centre_entry_index, prm_num_entries );
		return result.first ? ::std::make_optional( result.second ) : ::std::nullopt;
	}

} // namespace cath::scan::detail

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_STRIDE_CO_STRIDE_HPP
