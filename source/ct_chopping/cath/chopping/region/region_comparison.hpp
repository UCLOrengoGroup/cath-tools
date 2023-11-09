/// \file
/// \brief The region_comparison enum header

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

#ifndef CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_REGION_REGION_COMPARISON_HPP
#define CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_REGION_REGION_COMPARISON_HPP

namespace cath {

	/// \brief The different possible categories of result from comparing the locations of two regions
	///        that are assumed to annotate the same chain with the same indexing system
	enum class region_comparison : char {
		STRICTLY_BEFORE,      ///< The stop  of the first region comes strictly before the start of the second
		OVERLAPPINGLY_BEFORE, ///< The first region's start comes strictly first but its stop  comes somewhere between the second's first  and penultimate residues
		THE_SAME_AS,          ///< The first region's start and stop match those of the second region
		OVERLAPPINGLY_AFTER,  ///< The first region's stop  comes strictly last  but its start comes somewhere between the second's second and last        residues
		STRICTLY_AFTER,       ///< The start of the first region comes strictly after  the stop  of the second

		A_SUBSET_OF,          ///< The first region is a strict subset of the second

		A_SUPERSET_OF         ///< The first region is a strict superset of the second
	};

} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_REGION_REGION_COMPARISON_HPP
