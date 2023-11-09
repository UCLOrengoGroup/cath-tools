/// \file
/// \brief The constexpr_clamp header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_CONSTEXPR_CLAMP_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_CONSTEXPR_CLAMP_HPP

#include <stdexcept>

namespace cath::common {

	/// \brief Return the result of clamping prm_value within the range [ prm_low, prm_high ]
	template <typename T, typename U>
	constexpr T constexpr_clamp( const T &prm_value, ///< The value to be clamped
	                             const U &prm_low,   ///< The minimum value to which the prm_value must be clamped
	                             const U &prm_high   ///< The maximum value to which the prm_value must be clamped
	                             ) {
		return ( prm_low > prm_high ) ? throw ::std::logic_error("Unable to clamp to invalid range")
		                              : ( prm_value < prm_low  ) ? prm_low  :
		                                ( prm_value > prm_high ) ? prm_high :
		                                                           prm_value;
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_CONSTEXPR_CLAMP_HPP
