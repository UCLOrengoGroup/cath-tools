/// \file
/// \brief The res_pair_core functions header

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

#ifndef RES_PAIR_CORE_FUNCTIONS_H_INCLUDED
#define RES_PAIR_CORE_FUNCTIONS_H_INCLUDED

#include <boost/geometry/algorithms/comparable_distance.hpp>

#include "scan/detail/res_pair/res_pair_core.h"
#include "scan/detail/scan_type_aliases.h"
#include "structure/geometry/coord.h"
#include "structure/geometry/quat_rot.h"

namespace cath {
	namespace scan {

		namespace detail {

			/// \brief Calculate the distance between the quaternions for the two pairs' view_frames
			///        (the coordinate frames of their to_residues as seen from the coordinate frames of their from_residues)
			///
			/// This distance can be calculated a bit quicker than the angle and can be used in an exactly equivalent criterion
			///
			/// \relates res_pair_core
			inline frame_quat_rot_type distance_1_between_frames(const res_pair_core &arg_res_pair_a, ///< The first  res_pair_core
			                                                     const res_pair_core &arg_res_pair_b  ///< The second res_pair_core
			                                                     ) {
				const frame_quat_rot frame_a = arg_res_pair_a.get_frame();
				const frame_quat_rot frame_b = arg_res_pair_b.get_frame();
				return distance_1_between_quat_rots( frame_a, frame_b );
			}

			/// \brief Calculate the squared distance between the views of the two residue pairs
			///
			/// Each view is the location of the to_residue's carbon-beta atom as seen from
			/// the coordinate frame of the from_residue
			///
			/// \relates res_pair_core
			inline double squared_distance(const res_pair_core &arg_res_pair_a, ///< The first  res_pair_core
			                               const res_pair_core &arg_res_pair_b  ///< The second res_pair_core
			                               ) {
				return boost::geometry::comparable_distance(
					arg_res_pair_a.get_view(),
					arg_res_pair_b.get_view()
				);
			}

			/// \brief Get the (wrapped) difference between the two pairs' from_residue phi angles
			///
			/// \relates res_pair_core
			inline angle_type from_phi_angle_difference(const res_pair_core &arg_res_pair_a, ///< The first  res_pair_core
			                                            const res_pair_core &arg_res_pair_b  ///< The second res_pair_core
			                                            ) {
				return unshifted_wrapped_difference(
					arg_res_pair_a.get_from_phi_angle(),
					arg_res_pair_b.get_from_phi_angle()
				);
			}

			/// \brief Get the (wrapped) difference between the two pairs' from_residue psi angles
			///
			/// \relates res_pair_core
			inline angle_type from_psi_angle_difference(const res_pair_core &arg_res_pair_a, ///< The first  res_pair_core
			                                            const res_pair_core &arg_res_pair_b  ///< The second res_pair_core
			                                            ) {
				return unshifted_wrapped_difference(
					arg_res_pair_a.get_from_psi_angle(),
					arg_res_pair_b.get_from_psi_angle()
				);
			}

			/// \brief Get the (wrapped) difference between the two pairs' to_residue phi angles
			///
			/// \relates res_pair_core
			inline angle_type to_phi_angle_difference(const res_pair_core &arg_res_pair_a,   ///< The first  res_pair_core
			                                          const res_pair_core &arg_res_pair_b    ///< The second res_pair_core
			                                          ) {
				return unshifted_wrapped_difference(
					arg_res_pair_a.get_to_phi_angle(),
					arg_res_pair_b.get_to_phi_angle()
				);
			}

			/// \brief Get the (wrapped) difference between the two pairs' to_residue psi angles
			///
			/// \relates res_pair_core
			inline angle_type to_psi_angle_difference(const res_pair_core &arg_res_pair_a,   ///< The first  res_pair_core
			                                          const res_pair_core &arg_res_pair_b    ///< The second res_pair_core
			                                          ) {
				return unshifted_wrapped_difference(
					arg_res_pair_a.get_to_psi_angle(),
					arg_res_pair_b.get_to_psi_angle()
				);
			}

			/// \brief Get the maximum (wrapped) difference between the two pairs' from/to residue phi angles
			///
			/// \relates res_pair_core
			inline angle_type max_phi_angle_difference(const res_pair_core &arg_res_pair_a,   ///< The first  res_pair_core
			                                           const res_pair_core &arg_res_pair_b    ///< The second res_pair_core
			                                           ) {
				return std::max(
					from_phi_angle_difference( arg_res_pair_a, arg_res_pair_b ),
					  to_phi_angle_difference( arg_res_pair_a, arg_res_pair_b )
				);
			}

			/// \brief Get the maximum (wrapped) difference between the two pairs' from/to residue psi angles
			///
			/// \relates res_pair_core
			inline angle_type max_psi_angle_difference(const res_pair_core &arg_res_pair_a,   ///< The first  res_pair_core
			                                           const res_pair_core &arg_res_pair_b    ///< The second res_pair_core
			                                           ) {
				return std::max(
					from_psi_angle_difference( arg_res_pair_a, arg_res_pair_b ),
					  to_psi_angle_difference( arg_res_pair_a, arg_res_pair_b )
				);
			}

		}
	}
}

#endif
