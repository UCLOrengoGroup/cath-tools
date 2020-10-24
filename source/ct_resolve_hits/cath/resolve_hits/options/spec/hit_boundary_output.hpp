/// \file
/// \brief The hit_boundary_output header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_HIT_BOUNDARY_OUTPUT_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_HIT_BOUNDARY_OUTPUT_HPP

namespace cath {
	namespace rslv {

		/// \brief Whether the output of the hits' boundaries should have them trimmed
		enum class hit_boundary_output : bool {
			TRIMMED, ///< Output the trimmed version of the boundaries
			ORIG     ///< Output the original, untrimmed version of the boundaries
		};

		/// \brief Convert the specified bool output-trimmed-hits to the equivalent hit_boundary_output
		///
		/// \relates hit_boundary_output
		inline constexpr hit_boundary_output hit_boundary_output_of_output_trimmed_hits(const bool &prm_output_trimmed_hits ///< Whether the output of the hits' boundaries should have them trimmed
		                                                                                ) {
			return prm_output_trimmed_hits ? hit_boundary_output::TRIMMED
			                               : hit_boundary_output::ORIG;
		}

		/// \brief Convert the specified hit_boundary_output to a bool representing denoting
		///        whether the output of the hits' boundaries should have them trimmed
		///
		/// \relates hit_boundary_output
		inline constexpr bool means_output_trimmed_hits(const hit_boundary_output &prm_hit_boundary_output ///< The hit_boundary_output to query
		                                                ) {
			return ( prm_hit_boundary_output == hit_boundary_output::TRIMMED );
		}

	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_HIT_BOUNDARY_OUTPUT_HPP
