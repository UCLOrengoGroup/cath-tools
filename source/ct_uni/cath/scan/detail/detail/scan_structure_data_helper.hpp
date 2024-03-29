/// \file
/// \brief The scan_structure_data helper header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_DETAIL_SCAN_STRUCTURE_DATA_HELPER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_DETAIL_SCAN_STRUCTURE_DATA_HELPER_HPP

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "cath/scan/detail/res_pair/single_struc_res_pair_list.hpp"
#include "cath/scan/detail/stride/rep_strider.hpp"
#include "cath/scan/detail/stride/roled_scan_stride.hpp"

namespace cath::scan::detail::detail {

	/// \brief TODOCUMENT
	inline angle_type_vec make_scan_phi_angles(const protein &prm_protein ///< TODOCUMENT
	                                           ) {
		angle_type_vec results;
		results.reserve( prm_protein.get_length() );
		for (const auto &x : prm_protein) {
			results.emplace_back( geom::convert_angle_type<angle_base_type>( x.get_phi_angle() ) );
		}
		return results;
	}

	/// \brief TODOCUMENT
	inline angle_type_vec make_scan_psi_angles(const protein &prm_protein ///< TODOCUMENT
	                                           ) {
		angle_type_vec results;
		results.reserve( prm_protein.get_length() );
		for (const auto &x : prm_protein) {
			results.emplace_back( geom::convert_angle_type<angle_base_type>( x.get_psi_angle() ) );
		}
		return results;
	}

	/// \brief TODOCUMENT
	inline view_type_vec make_scan_view_coords(const protein &prm_protein ///< TODOCUMENT
	                                           ) {
		view_type_vec results;
		results.reserve( prm_protein.get_length() );
		for (const auto &x : prm_protein) {
			results.emplace_back( x.get_carbon_beta_coord() );
		}
		return results;
	}

	/// \brief TODOCUMENT
	inline frame_quat_rot_vec make_scan_frame_quat_rots(const protein &prm_protein ///< TODOCUMENT
	                                                    ) {
		frame_quat_rot_vec results;
		results.reserve( prm_protein.get_length() );
		for (const auto &x : prm_protein) {
			results.emplace_back( geom::make_quat_rot_from_rotation<frame_quat_rot_type>( x.get_frame() ) );
		}
		return results;
	}

	/// \brief TODOCUMENT
	inline single_struc_res_pair_vec build_single_rep_pairs(const protein &prm_protein ///< TODOCUMENT
	                                                        ) {
		const auto num_residues     = debug_unwarned_numeric_cast<index_type>( prm_protein.get_length() );
		const auto scan_phi_angles  = make_scan_phi_angles     ( prm_protein );
		const auto scan_psi_angles  = make_scan_psi_angles     ( prm_protein );
		const auto scan_view_coords = make_scan_view_coords    ( prm_protein );
		const auto scan_frames      = make_scan_frame_quat_rots( prm_protein );
		single_struc_res_pair_vec results;
		results.reserve( num_residues * num_residues );
		const auto all_residues_range   = boost::irange<index_type>( 0, num_residues );
		for (const auto &x : common::cross( all_residues_range, all_residues_range ) ) {
			const auto &from_res_index = std::get<0>( x );
			const auto &to_res_index   = std::get<1>( x );
			results.emplace_back(
				rotate_copy(
					scan_frames     [ from_res_index ],
					scan_view_coords[ to_res_index   ] - scan_view_coords[ from_res_index ]
				),
				rotation_between_rotations(
					scan_frames[ from_res_index ],
					scan_frames[ to_res_index   ]
				),
				scan_phi_angles[ from_res_index ],
				scan_psi_angles[ from_res_index ],
				scan_phi_angles[ to_res_index   ],
				scan_psi_angles[ to_res_index   ],
				from_res_index,
				to_res_index
			);
		}
		return results;
	}

	/// \brief TODOCUMENT
	inline void add_stride_neighbours(single_struc_res_pair_list      &prm_neighbours,       ///< TODOCUMENT
	                                  const index_type                &prm_num_residues,     ///< TODOCUMENT
	                                  const single_struc_res_pair_vec &all_single_res_pairs, ///< TODOCUMENT
	                                  const index_type                &prm_from_index,       ///< TODOCUMENT
	                                  const index_type                &prm_to_index,         ///< TODOCUMENT
	                                  const index_type                &prm_from_co_stride,   ///< TODOCUMENT
	                                  const index_type                &prm_to_co_stride      ///< TODOCUMENT
	                                  ) {
		const auto from_stride_size  = num_in_stride_neighbour_range( prm_from_co_stride );
		const auto to_stride_size    = num_in_stride_neighbour_range( prm_to_co_stride   );
		const auto from_stride_range = boost::irange<index_type>( 0, from_stride_size );
		const auto to_stride_range   = boost::irange<index_type>( 0, to_stride_size   );

		prm_neighbours.reserve( from_stride_size * to_stride_size );

		for(const auto &x : common::cross( from_stride_range, to_stride_range )) {
			const auto from = entry_index_of_stride_neighbour_index( std::get<0>( x ), prm_from_co_stride, prm_from_index, prm_num_residues );
			const auto to   = entry_index_of_stride_neighbour_index( std::get<1>( x ), prm_to_co_stride,   prm_to_index,   prm_num_residues );
			if ( from && to && ( *from != *to ) ) {
				prm_neighbours.push_back( all_single_res_pairs[ ( (*from) * prm_num_residues ) + (*to) ] );
			}
			else {
				prm_neighbours.emplace_back();
			}
		};
	}

	/// \brief TODOCUMENT
	inline single_struc_res_pair_list_vec build_rep_sets(const protein           &prm_protein,          ///< TODOCUMENT
	                                                     const roled_scan_stride &prm_roled_scan_stride ///< TODOCUMENT
	                                                     ) {
		const auto  num_residues         = debug_unwarned_numeric_cast<index_type>( prm_protein.get_length() );
		const auto &all_single_res_pairs = build_single_rep_pairs( prm_protein );

		const auto from_range = get_indices_range( get_this_from_strider( prm_roled_scan_stride ), num_residues );
		const auto to_range   = get_indices_range( get_this_to_strider  ( prm_roled_scan_stride ), num_residues );

		single_struc_res_pair_list_vec results;
		results.reserve( from_range.size() * to_range.size() );

		for (const auto &x : common::cross( from_range, to_range ) ) {
//			std::cerr << "Building rep: [" << std::get<0>( x ) << ", " << std::get<1>( x ) << "]\n";
			results.emplace_back();
			detail::add_stride_neighbours(
				results.back(),
				num_residues,
				all_single_res_pairs,
				std::get<0>( x ),
				std::get<1>( x ),
				from_co_stride( prm_roled_scan_stride ),
				to_co_stride  ( prm_roled_scan_stride )
			);
		}
		return results;
	}

} // namespace cath::scan::detail::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_DETAIL_SCAN_STRUCTURE_DATA_HELPER_HPP
