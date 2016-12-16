/// \file
/// \brief The dssp_skip_policy class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_DSSP_SKIP_POLICY_DSSP_SKIP_POLICY_H
#define _CATH_TOOLS_SOURCE_FILE_DSSP_SKIP_POLICY_DSSP_SKIP_POLICY_H

namespace cath {
	namespace file {

		/// \brief Whether phi/psi angles should be skipped for residues that DSSP would skip
		///        (and for the preceding residue's psi and following residues phi)
		///
		/// Note: this doesn't affect objective breaks in the chain (where there is a big gap
		/// between one residue's C and the next residues N or where the chain label switches)
		/// for which the phi/psi angles are always skipped
		enum class dssp_skip_angle_skipping : bool {
			BREAK_ANGLES,     ///< Break the chain of phi/psi angles at residues DSSP would skip
			DONT_BREAK_ANGLES ///< Process phi/psi angles normally for residues that that DSSP would skip
			                  ///< (but still break phi/psi angles at genuine breaks in the chain)
		};

		/// \brief Whether to skip residues that DSSP would skip
		enum class dssp_skip_res_skipping : bool {
			SKIP,     ///< Skip residues that DSSP would skip
			DONT_SKIP ///< Don't skip residues that DSSP would skip
		};

		/// \brief Policy for handling residues that DSSP would skip, combining dssp_skip_angle_skipping and dssp_skip_res_skipping
		///
		/// Note that there is no option to skip but not break angles. This is completely impossible
		/// but doesn't fit in well with the current code operation so should only be added if needed.
		enum class dssp_skip_policy : char {
			SKIP__BREAK_ANGLES,          ///< Skip residues that DSSP would skip and break connecting phi/psi angles of neighbouring residues
			DONT_SKIP__BREAK_ANGLES,     ///< Don't skip residues that DSSP would skip but do break phi/psi angles of those residues and connecting phi/psi angles of neighbouring residues
			DONT_SKIP__DONT_BREAK_ANGLES ///< Don't skip residues that DSSP would skip and don't break phi/psi angles
		};

		/// \brief Extract the angle-skipping aspect of the specified dssp_skip_policy
		constexpr dssp_skip_angle_skipping angle_skipping_of_dssp_skip_policy(const dssp_skip_policy &arg_dssp_skip_policy ///< The dssp_skip_policy to query
		                                                                      ) {
			return ( arg_dssp_skip_policy == dssp_skip_policy::DONT_SKIP__DONT_BREAK_ANGLES )
				? dssp_skip_angle_skipping::DONT_BREAK_ANGLES
				: dssp_skip_angle_skipping::BREAK_ANGLES;
		}

		/// \brief Extract the residue-skipping aspect of the specified dssp_skip_policy
		constexpr dssp_skip_res_skipping res_skipping_of_dssp_skip_policy(const dssp_skip_policy &arg_dssp_skip_policy ///< The dssp_skip_policy to query
		                                                                  ) {
			return ( arg_dssp_skip_policy == dssp_skip_policy::SKIP__BREAK_ANGLES )
				? dssp_skip_res_skipping::SKIP
				: dssp_skip_res_skipping::DONT_SKIP;
		}

	} // namespace file
} // namespace cath

#endif
