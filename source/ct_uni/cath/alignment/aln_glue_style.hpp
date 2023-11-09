/// \file
/// \brief The aln_glue_style header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALN_GLUE_STYLE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALN_GLUE_STYLE_HPP

namespace cath::align {

	/// \brief The style of approach to glueing alignments together
	enum class aln_glue_style : char {
		INCREMENTALLY_WITH_PAIR_REFINING, ///< Incrementally glue in alignments that add one entry and then refine just that new alignment
		SIMPLY,                           ///< Just glue the alignments together
		WITH_HEAVY_REFINING               ///< Do oodles of really heavy alignment
	};

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALN_GLUE_STYLE_HPP
