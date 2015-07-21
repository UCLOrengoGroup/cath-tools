/// \file
/// \brief The quad_criteria_result class header

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

#ifndef QUAD_CRITERIA_RESULT_H_INCLUDED
#define QUAD_CRITERIA_RESULT_H_INCLUDED

#include <iosfwd>

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief TODOCUMENT
			///
			/// This includes whether the quad in quest
			enum class quad_criteria_result : size_t {
				PASS,                      ///< TODOCUMENT
				QUERY_FAILS_SINGLE_CHECKS, ///< TODOCUMENT
				INDEX_FAILS_SINGLE_CHECKS, ///< TODOCUMENT
				FAILS_VIEW_CHECK,          ///< TODOCUMENT
				FAILS_PHI_CHECK,           ///< TODOCUMENT
				FAILS_PSI_CHECK,           ///< TODOCUMENT
				FAILS_FRAME_CHECK,         ///< TODOCUMENT
				FAILS_QUAD_CHECKS,         ///< TODOCUMENT
				HAS_NO_REP                 ///< TODOCUMENT
			};

			std::string to_string(const quad_criteria_result &);
			std::ostream & operator<<(std::ostream &,
			                          const quad_criteria_result &);

		}
	}
}

#endif
