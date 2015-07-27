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

#ifndef CONSTEXPR_CLAMP_H_INCLUDED
#define CONSTEXPR_CLAMP_H_INCLUDED

namespace cath {
	namespace common {

		/// \brief Return the result of clamping arg_value within the range [ arg_low, arg_high ]
		template <typename T, typename U>
		constexpr const T constexpr_clamp(const T &arg_value, ///< The value to be clamped
		                                  const U &arg_low,   ///< The minimum value to which the arg_value must be clamped
		                                  const U &arg_high   ///< The maximum value to which the arg_value must be clamped
		                                  ) {
			return ( arg_low > arg_high ) ? throw std::logic_error("Unable to clamp to invalid range")
			                              : ( arg_value < arg_low  ) ? arg_low  :
			                                ( arg_value > arg_high ) ? arg_high :
			                                                           arg_value;
		}

	}
}

#endif
