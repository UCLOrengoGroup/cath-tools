/// \file
/// \brief compare_upper_cell_result enum

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 1989, Orengo Group, University College London
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

#ifndef _CATH_TOOLS_SOURCE_SSAP_COMPARE_UPPER_CELL_RESULT_H
#define _CATH_TOOLS_SOURCE_SSAP_COMPARE_UPPER_CELL_RESULT_H

namespace cath {

	/// \brief The result of a compare_upper_cell() call
	enum class compare_upper_cell_result : char {
		ZERO,                     ///< The lower alignment scored zero and so won't be added to the upper matrix
		NON_ZERO_BELOW_THRESHOLD, ///< The lower alignment scored more than zero but not enough to meet the MIN_LOWER_MAT_RES_SCORE threshold,
		                          ///< so it won't be added to the upper matrix
		SCORED                    ///< The lower alignment scored enough to meet the MIN_LOWER_MAT_RES_SCORE threshold and
		                          ///< will be added to the upper matrix
	};
} // namespace cath

#endif
