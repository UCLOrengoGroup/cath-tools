/// \file
/// \brief The res_pair_core functions header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_FUNCTIONS_RES_PAIR_CORE_FUNCTIONS_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_FUNCTIONS_RES_PAIR_CORE_FUNCTIONS_HPP

#include <boost/geometry/algorithms/comparable_distance.hpp>

#include "cath/scan/detail/res_pair/res_pair_core.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/quat_rot.hpp"

namespace cath::scan::detail {

	/// \brief Calculate the distance between the quaternions for the two pairs' view_frames
	///        (the coordinate frames of their to_residues as seen from the coordinate frames of their from_residues)
	///
	/// This distance can be calculated a bit quicker than the angle and can be used in an exactly equivalent criterion
	///
	/// \relates res_pair_core
	inline frame_quat_rot_type distance_1_between_frames(const res_pair_core &prm_res_pair_a, ///< The first  res_pair_core
	                                                     const res_pair_core &prm_res_pair_b  ///< The second res_pair_core
	                                                     ) {
		return distance_1_between_quat_rots( prm_res_pair_a.get_frame(), prm_res_pair_b.get_frame() );
	}

	/// \brief Calculate the squared distance between the views of the two residue pairs
	///
	/// Each view is the location of the to_residue's carbon-beta atom as seen from
	/// the coordinate frame of the from_residue
	///
	/// \relates res_pair_core
	inline double squared_distance(const res_pair_core &prm_res_pair_a, ///< The first  res_pair_core
	                               const res_pair_core &prm_res_pair_b  ///< The second res_pair_core
	                               ) {
		return boost::geometry::comparable_distance(
			prm_res_pair_a.get_view(),
			prm_res_pair_b.get_view()
		);
	}

	/// \brief Get the (wrapped) difference between the two pairs' from_residue phi angles
	///
	/// \relates res_pair_core
	inline angle_type from_phi_angle_difference(const res_pair_core &prm_res_pair_a, ///< The first  res_pair_core
	                                            const res_pair_core &prm_res_pair_b  ///< The second res_pair_core
	                                            ) {
		return unshifted_wrapped_difference(
			prm_res_pair_a.get_from_phi_angle(),
			prm_res_pair_b.get_from_phi_angle()
		);
	}

	/// \brief Get the (wrapped) difference between the two pairs' from_residue psi angles
	///
	/// \relates res_pair_core
	inline angle_type from_psi_angle_difference(const res_pair_core &prm_res_pair_a, ///< The first  res_pair_core
	                                            const res_pair_core &prm_res_pair_b  ///< The second res_pair_core
	                                            ) {
		return unshifted_wrapped_difference(
			prm_res_pair_a.get_from_psi_angle(),
			prm_res_pair_b.get_from_psi_angle()
		);
	}

	/// \brief Get the (wrapped) difference between the two pairs' to_residue phi angles
	///
	/// \relates res_pair_core
	inline angle_type to_phi_angle_difference(const res_pair_core &prm_res_pair_a,   ///< The first  res_pair_core
	                                          const res_pair_core &prm_res_pair_b    ///< The second res_pair_core
	                                          ) {
		return unshifted_wrapped_difference(
			prm_res_pair_a.get_to_phi_angle(),
			prm_res_pair_b.get_to_phi_angle()
		);
	}

	/// \brief Get the (wrapped) difference between the two pairs' to_residue psi angles
	///
	/// \relates res_pair_core
	inline angle_type to_psi_angle_difference(const res_pair_core &prm_res_pair_a,   ///< The first  res_pair_core
	                                          const res_pair_core &prm_res_pair_b    ///< The second res_pair_core
	                                          ) {
		return unshifted_wrapped_difference(
			prm_res_pair_a.get_to_psi_angle(),
			prm_res_pair_b.get_to_psi_angle()
		);
	}

	/// \brief Get the maximum (wrapped) difference between the two pairs' from/to residue phi angles
	///
	/// \relates res_pair_core
	inline angle_type max_phi_angle_difference(const res_pair_core &prm_res_pair_a,   ///< The first  res_pair_core
	                                           const res_pair_core &prm_res_pair_b    ///< The second res_pair_core
	                                           ) {
		return std::max(
			from_phi_angle_difference( prm_res_pair_a, prm_res_pair_b ),
			  to_phi_angle_difference( prm_res_pair_a, prm_res_pair_b )
		);
	}

	/// \brief Get the maximum (wrapped) difference between the two pairs' from/to residue psi angles
	///
	/// \relates res_pair_core
	inline angle_type max_psi_angle_difference(const res_pair_core &prm_res_pair_a,   ///< The first  res_pair_core
	                                           const res_pair_core &prm_res_pair_b    ///< The second res_pair_core
	                                           ) {
		return std::max(
			from_psi_angle_difference( prm_res_pair_a, prm_res_pair_b ),
			  to_psi_angle_difference( prm_res_pair_a, prm_res_pair_b )
		);
	}

} // namespace cath::scan::detail

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_FUNCTIONS_RES_PAIR_CORE_FUNCTIONS_HPP
