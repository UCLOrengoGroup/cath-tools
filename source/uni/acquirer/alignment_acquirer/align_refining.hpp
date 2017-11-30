/// \file
/// \brief The align_refining class header

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

#ifndef _CATH_TOOLS_SOURCE_ACQUIRER_ALIGNMENT_ACQUIRER_ALIGN_REFINING_H
#define _CATH_TOOLS_SOURCE_ACQUIRER_ALIGNMENT_ACQUIRER_ALIGN_REFINING_H

namespace cath {
	namespace align {

		/// \brief Represent how much refining should be done to an alignment
		///        (typically as part of the alignment being acquired or immediately afterwards)
		enum class align_refining {
			NO,    ///< No refining should be performed
			LIGHT, ///< At most, light refining should be performed
			HEAVY  ///< Heavy, slow, expensive refining should be performed
		};

	} // namespace align
} // namespace cath

#endif
