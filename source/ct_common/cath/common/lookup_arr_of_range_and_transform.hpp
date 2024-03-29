/// \file
/// \brief The lookup_arr_of_range_and_transform header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_LOOKUP_ARR_OF_RANGE_AND_TRANSFORM_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_LOOKUP_ARR_OF_RANGE_AND_TRANSFORM_HPP

#include <array>
#include <tuple>
#include <utility>

#include "cath/common/cpp20/constexpr_invoke.hpp"

namespace cath::common {

	template <typename Rng, typename Fn>
	constexpr auto lookup_arr_of_range_and_transform( Rng &&prm_rng, Fn &&prm_fn ) {
		return ::std::apply(
		  [ & ]( auto &&...prms ) {
			  return ::std::array{ ::std::pair{ prms, constexpr_invoke( prm_fn, prms ) }... };
		  },
		  prm_rng );
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_LOOKUP_ARR_OF_RANGE_AND_TRANSFORM_HPP
