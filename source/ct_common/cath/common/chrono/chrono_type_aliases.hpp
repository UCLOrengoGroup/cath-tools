/// \file
/// \brief The chrono type_aliases header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CHRONO_CHRONO_TYPE_ALIASES_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CHRONO_CHRONO_TYPE_ALIASES_HPP

#include <chrono>
#include <vector>

namespace cath {
	/// \brief TODOCUMENT
	using hrc_time_point = std::chrono::high_resolution_clock::time_point;

	/// \brief TODOCUMENT
	using hrc_time_point_vec = std::vector<hrc_time_point>;

	/// \brief TODOCUMENT
	using hrc_duration = std::chrono::high_resolution_clock::duration;

	/// \brief TODOCUMENT
	using hrc_duration_vec = std::vector<hrc_duration>;
} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CHRONO_CHRONO_TYPE_ALIASES_HPP
