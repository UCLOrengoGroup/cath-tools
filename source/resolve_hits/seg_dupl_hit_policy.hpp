/// \file
/// \brief The seg_dupl_hit_policy class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_SEG_DUPL_HIT_POLICY_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_SEG_DUPL_HIT_POLICY_H

namespace cath {
	namespace rslv {

		/// \brief What to do with pairs of hits with identical segments
		enum class seg_dupl_hit_policy : bool {
			PRESERVE, ///< Preserve both hits
			PRUNE     ///< Prune the worse (non-better) of the two hits
		};

	} // namespace rslv
} // namespace cath

#endif
