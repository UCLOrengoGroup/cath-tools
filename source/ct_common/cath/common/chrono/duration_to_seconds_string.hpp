/// \file
/// \brief The duration_to_seconds_string header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CHRONO_DURATION_TO_SECONDS_STRING_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CHRONO_DURATION_TO_SECONDS_STRING_HPP

#include <chrono>
#include <string>

namespace cath::common {

	/// \brief TODOCUMENT
	template <typename DURN>
	inline double durn_to_seconds_double(const DURN &prm_duration ///< The duration to convert
	                                     ) {
		return std::chrono::duration_cast<std::chrono::duration<double>>(
			prm_duration
		).count();
	}

	/// \brief TODOCUMENT
	template <typename DURN>
	inline std::string durn_to_seconds_string(const DURN &prm_duration ///< The duration to convert
	                                          ) {
		return std::to_string(
			std::chrono::duration_cast<std::chrono::duration<double>>(
				prm_duration
			).count()
		) + " seconds";
	}

	/// \brief TODOCUMENT
	template <typename DURN>
	inline double durn_to_rate_per_second(const DURN &prm_duration ///< The duration to convert
	                                      ) {
		return 1.0 / durn_to_seconds_double( prm_duration );
	}

	/// \brief TODOCUMENT
	template <typename DURN>
	inline std::string durn_to_rate_per_second_string(const DURN &prm_duration ///< The duration to convert
	                                                  ) {
		return std::to_string( durn_to_rate_per_second( prm_duration ) ) + " / second";
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CHRONO_DURATION_TO_SECONDS_STRING_HPP

