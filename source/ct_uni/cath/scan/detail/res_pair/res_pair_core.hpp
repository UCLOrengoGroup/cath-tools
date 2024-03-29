/// \file
/// \brief The res_pair_core class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_RES_PAIR_CORE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_RES_PAIR_CORE_HPP

#include "cath/scan/detail/res_pair/functions/res_index_pair_functions.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"
#include "cath/ssap/context_res.hpp"
#include "cath/structure/geometry/angle.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/quat_rot.hpp"
#include "cath/structure/protein/residue.hpp"

#include <iosfwd>

namespace cath::scan::detail {

	/// \brief Describe the core info for a pair of from / to residues for quick searching
	///
	/// This is used in the implementation of single_struc_res_pair and multi_struc_res_rep_pair
	///
	/// Invariants:
	///  * All four from/to phi/psi angles must be shifted to the same (standard) angle range
	class res_pair_core final {
	private:
		/// \brief The view of the to_residue from the from_residue
		/// (ie the vector from the from_residue to the to_residue, in the coordinate frame of the from_residue)
		view_type view;

		/// \brief The coordinate frame of the from_residue, as determined by its core atoms
		frame_quat_rot frame;

		/// \brief The phi angle of the from_residue
		angle_type from_phi_angle;

		/// \brief The psi angle of the from_residue
		angle_type from_psi_angle;

		/// \brief The phi angle of the to_residue
		angle_type to_phi_angle;

		/// \brief The psi angle of the to_residue
		angle_type to_psi_angle;

	public:
		res_pair_core();
		res_pair_core(view_type,
		              frame_quat_rot,
		              angle_type,
		              angle_type,
		              angle_type,
		              angle_type);
		[[nodiscard]] const view_type &     get_view() const;
		[[nodiscard]] const frame_quat_rot &get_frame() const;
		[[nodiscard]] const angle_type &    get_from_phi_angle() const;
		[[nodiscard]] const angle_type &    get_from_psi_angle() const;
		[[nodiscard]] const angle_type &    get_to_phi_angle() const;
		[[nodiscard]] const angle_type &    get_to_psi_angle() const;
	};

	/// \brief Ctor to create dummy res_pair
	inline res_pair_core::res_pair_core() : frame( geom::make_identity_quat_rot<frame_quat_rot_type>() ),
	                                        from_phi_angle( 0.0 ),
	                                        from_psi_angle( 0.0 ),
	                                        to_phi_angle  ( 0.0 ),
	                                        to_psi_angle  ( 0.0 ) {
	}

	/// \brief Ctor to populate all members
	inline res_pair_core::res_pair_core(view_type      prm_view,           ///< The view of the to_residue from the from_residue
	                                    frame_quat_rot prm_frame,          ///< The coordinate frame of the from_residue, as determined by its core atoms
	                                    angle_type     prm_from_phi_angle, ///< The phi angle of the from_residue
	                                    angle_type     prm_from_psi_angle, ///< The psi angle of the from_residue
	                                    angle_type     prm_to_phi_angle,   ///< The phi angle of the to_residue
	                                    angle_type     prm_to_psi_angle    ///< The psi angle of the to_residue
	                                    ) : view           { std::move( prm_view           ) },
	                                        frame          { std::move( prm_frame          ) },
	                                        from_phi_angle { std::move( prm_from_phi_angle ) },
	                                        from_psi_angle { std::move( prm_from_psi_angle ) },
	                                        to_phi_angle   { std::move( prm_to_phi_angle   ) },
	                                        to_psi_angle   { std::move( prm_to_psi_angle   ) } {
	}

	/// \brief Getter for view
	inline const view_type & res_pair_core::get_view() const {
		return view;
	}

	/// \brief Getter for frame
	inline const frame_quat_rot & res_pair_core::get_frame() const {
		return frame;
	}

	/// \brief Getter for from_phi_angle
	inline const angle_type & res_pair_core::get_from_phi_angle() const {
		return from_phi_angle;
	}

	/// \brief Getter for from_psi_angle
	inline const angle_type & res_pair_core::get_from_psi_angle() const {
		return from_psi_angle;
	}

	/// \brief Getter for to_phi_angle
	inline const angle_type & res_pair_core::get_to_phi_angle() const {
		return to_phi_angle;
	}

	/// \brief Getter for to_psi_angle
	inline const angle_type & res_pair_core::get_to_psi_angle() const {
		return to_psi_angle;
	}

	/// \brief Convenience function to get the x component of the view in the specified res_pair_core
	///
	/// \relates res_pair_core
	inline const view_base_type & get_view_x(const res_pair_core &prm_res_pair ///< The res_pair_core to query
	                                          ) {
	        return prm_res_pair.get_view().get<0>();
	}

	/// \brief Convenience function to get the y component of the view in the specified res_pair_core
	///
	/// \relates res_pair_core
	inline const view_base_type & get_view_y(const res_pair_core &prm_res_pair ///< The res_pair_core to query
	                                          ) {
	        return prm_res_pair.get_view().get<1>();
	}

	/// \brief Convenience function to get the z component of the view in the specified res_pair_core
	///
	/// \relates res_pair_core
	inline const view_base_type & get_view_z(const res_pair_core &prm_res_pair ///< The res_pair_core to query
	                                          ) {
	        return prm_res_pair.get_view().get<2>();
	}

	/// \brief TODOCUMENT
	///
	/// \relates res_pair_core
	inline res_pair_core make_res_pair_core(const residue &prm_from_residue, ///< TODOCUMENT
	                                        const residue &prm_to_residue    ///< TODOCUMENT
	                                        ) {
		return {
			view_type{ view_vector_of_residue_pair( prm_from_residue, prm_to_residue ) },
			geom::make_quat_rot_from_rotation<frame_quat_rot_type>(
				view_frame( prm_from_residue, prm_to_residue )
			),
			geom::convert_angle_type<angle_base_type>( prm_from_residue.get_phi_angle() ).quick_shift(),
			geom::convert_angle_type<angle_base_type>( prm_from_residue.get_psi_angle() ).quick_shift(),
			geom::convert_angle_type<angle_base_type>( prm_to_residue.get_phi_angle()   ).quick_shift(),
			geom::convert_angle_type<angle_base_type>( prm_to_residue.get_psi_angle()   ).quick_shift()
		};
	}

	std::ostream & operator<<(std::ostream &,
	                          const res_pair_core &);

} // namespace cath::scan::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_RES_PAIR_CORE_HPP
