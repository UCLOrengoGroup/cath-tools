/// \file
/// \brief The array_of_first_n header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_MAKE_TYPE_OF_FIRST_N_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_MAKE_TYPE_OF_FIRST_N_HPP

#include <array>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"

namespace cath::common {

	namespace detail {

		/// Consider making this move values out of an rvalue range
		template <typename T, size_t N, typename Rng, size_t... Index>
		constexpr T make_type_of_first_n_impl( Rng &&prm_range, ::std::index_sequence<Index...> ) {
			return T{ prm_range[ Index ]... };
		}

	} // namespace detail

	/// \brief Return an array of size N with the first N elements of the specified range
	///
	/// \param prm_range The range from which the array should be built
	template <typename T, size_t N, typename Rng>
	constexpr T make_type_of_first_n( Rng &&prm_range ) {
		return detail::make_type_of_first_n_impl<T, N>( ::std::forward<Rng>( prm_range ), ::std::make_index_sequence<N>{} );
	}

	/// \brief Return an array of size N with the first N elements of the specified range
	///
	/// \param prm_range The range from which the array should be built
	template <size_t N, typename Rng>
	constexpr auto make_array_of_first_n( Rng &&prm_range ) {
		return detail::make_type_of_first_n_impl<::std::array<::cath::common::range_value_t<Rng>, N>, N>(
		  ::std::forward<Rng>( prm_range ), ::std::make_index_sequence<N>{} );
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_MAKE_TYPE_OF_FIRST_N_HPP
