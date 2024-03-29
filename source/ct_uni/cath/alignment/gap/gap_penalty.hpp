/// \file
/// \brief The gap_penalty class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_GAP_GAP_PENALTY_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_GAP_GAP_PENALTY_HPP

#include "cath/common/type_aliases.hpp"

namespace cath::align::gap {

	/// \brief TODOCUMENT
	///
	/// This is an affine gap penalty (ie permits a penalty for opening and another for extending)
	class gap_penalty final {
	private:
		/// \brief TODOCUMENT
		score_type open_gap_penalty;

		/// \brief TODOCUMENT
		score_type extend_gap_penalty;

	public:
		gap_penalty(const score_type &,
		            const score_type &);

		[[nodiscard]] score_type get_open_gap_penalty() const;
		[[nodiscard]] score_type get_extend_gap_penalty() const;
	};

} // namespace cath::align::gap

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_GAP_GAP_PENALTY_HPP
