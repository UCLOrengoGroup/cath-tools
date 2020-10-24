/// \file
/// \brief The view_cache_index_entry protein functions header.

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
///
/// This contains functions that are equivalent to the many of the non-member view_cache_index_entry functions
/// but taking protein and size_size_pair arguments rather than view_cache_index_entry arguments.
///
/// This is helpful in implementing equivalent functionality with and without view_cache_index code
/// so that the results can be compared etc.

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_VCIE_PROTEIN_FUNCTIONS_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_VCIE_PROTEIN_FUNCTIONS_HPP

#include "cath/ssap/context_res.hpp"
#include "cath/structure/geometry/rotation.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/view_cache/index/detail/view_cache_index_type_aliases.hpp"

#include <utility>

namespace cath { class protein; }
namespace cath { namespace geom { class rotation; } }

namespace cath {
	namespace index {
		namespace detail {

			/// \brief Calculate the view_frame of the pair
			///        (ie the coordinate frame of the to_residue's atoms in terms of the coordinate from of the from_residue's atoms)
			inline cath::geom::rotation view_frame(const protein &prm_protein,    ///< The protein containing the two residues
			                                       const size_t  &prm_from_index, ///< The index of the from_residue
			                                       const size_t  &prm_to_index    ///< The index of the to_residue
			                                       ) {
				const residue              &from_a_res   = prm_protein.get_residue_ref_of_index( prm_from_index );
				const residue              &to_a_res     = prm_protein.get_residue_ref_of_index( prm_to_index   );
				const cath::geom::rotation &from_a_frame = from_a_res.get_frame();
				const cath::geom::rotation &to_a_frame   = to_a_res.get_frame();
				return rotation_between_rotations( from_a_frame, to_a_frame );
			}

			/// \brief Calculate the squared distance between the views of the two residue pairs
			///
			/// Each view is the location of the to_residue's carbon-beta atom as seen from
			/// the coordinate frame of the from_residue
			inline double squared_distance(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                               const size_size_pair &prm_indices_b, ///< The indices of the from/to residues in the second protein
			                               const protein        &prm_protein_a, ///< The first  protein
			                               const protein        &prm_protein_b  ///< The second protein
			                               ) {
				const residue    &from_a_res = prm_protein_a.get_residue_ref_of_index( prm_indices_a.first  );
				const residue    &to_a_res   = prm_protein_a.get_residue_ref_of_index( prm_indices_a.second );
				const residue    &from_b_res = prm_protein_b.get_residue_ref_of_index( prm_indices_b.first  );
				const residue    &to_b_res   = prm_protein_b.get_residue_ref_of_index( prm_indices_b.second );

				const geom::coord view_a     = view_vector_of_residue_pair( from_a_res, to_a_res );
				const geom::coord view_b     = view_vector_of_residue_pair( from_b_res, to_b_res );

				return squared_distance_between_points( view_a, view_b );
			}

			/// \brief Return whether the two pairs of residues are both in the same direction
			///        (ie both have from_index < to_index or both have from_index > to_index)
			inline bool same_direction(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                           const size_size_pair &prm_indices_b  ///< The indices of the from/to residues in the second protein
			                           ) {
				const bool a_increases = ( prm_indices_a.first < prm_indices_a.second );
				const bool b_increases = ( prm_indices_b.first < prm_indices_b.second );
				return ( a_increases == b_increases );
			}

			/// \brief Calculate the minimum of the two pairs' absolute differences between their from_index and their to_index
			inline size_t min_index_difference(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                                   const size_size_pair &prm_indices_b  ///< The indices of the from/to residues in the second protein
			                                   ) {
				return std::min(
					common::difference( prm_indices_a.first, prm_indices_a.second ),
					common::difference( prm_indices_b.first, prm_indices_b.second )
				);
			}

			/// \brief Calculate the angle between the two pairs' view_frames
			///        (the coordinate frames of their to_residues as seen from the coordinate frames of their from_residues)
			inline angle_type angle_between_frames(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                                       const size_size_pair &prm_indices_b, ///< The indices of the from/to residues in the second protein
			                                       const protein        &prm_protein_a, ///< The first  protein
			                                       const protein        &prm_protein_b  ///< The second protein
			                                       ) {
				const geom::rotation a_frame = view_frame( prm_protein_a, prm_indices_a.first, prm_indices_a.second );
				const geom::rotation b_frame = view_frame( prm_protein_b, prm_indices_b.first, prm_indices_b.second );
				return geom::convert_angle_type<angle_base_type>( angle_between_rotations( a_frame, b_frame ) );
			}

			/// \brief Get the (wrapped) difference between the two pairs' from_residue phi angles
			inline angle_type from_phi_angle_difference(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                                            const size_size_pair &prm_indices_b, ///< The indices of the from/to residues in the second protein
			                                            const protein        &prm_protein_a, ///< The first  protein
			                                            const protein        &prm_protein_b  ///< The second protein
			                                            ) {
				const residue    &from_a_res = prm_protein_a.get_residue_ref_of_index( prm_indices_a.first );
				const residue    &from_b_res = prm_protein_b.get_residue_ref_of_index( prm_indices_b.first );
				const angle_type &from_phi_a = geom::convert_angle_type<angle_base_type>( from_a_res.get_phi_angle() );
				const angle_type &from_phi_b = geom::convert_angle_type<angle_base_type>( from_b_res.get_phi_angle() );
				return wrapped_difference( from_phi_a, from_phi_b );
			}

			/// \brief Get the (wrapped) difference between the two pairs' from_residue psi angles
			inline angle_type from_psi_angle_difference(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                                            const size_size_pair &prm_indices_b, ///< The indices of the from/to residues in the second protein
			                                            const protein        &prm_protein_a, ///< The first  protein
			                                            const protein        &prm_protein_b  ///< The second protein
			                                            ) {
				const residue &from_a_res = prm_protein_a.get_residue_ref_of_index( prm_indices_a.first );
				const residue &from_b_res = prm_protein_b.get_residue_ref_of_index( prm_indices_b.first );
				const auto     from_psi_a = geom::convert_angle_type<angle_base_type>( from_a_res.get_psi_angle() );
				const auto     from_psi_b = geom::convert_angle_type<angle_base_type>( from_b_res.get_psi_angle() );
				return wrapped_difference( from_psi_a, from_psi_b );
			}

			/// \brief Get the (wrapped) difference between the two pairs' to_residue phi angles
			inline angle_type to_phi_angle_difference(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                                          const size_size_pair &prm_indices_b, ///< The indices of the from/to residues in the second protein
			                                          const protein        &prm_protein_a, ///< The first  protein
			                                          const protein        &prm_protein_b  ///< The second protein
			                                          ) {
				const residue &to_a_res = prm_protein_a.get_residue_ref_of_index( prm_indices_a.second );
				const residue &to_b_res = prm_protein_b.get_residue_ref_of_index( prm_indices_b.second );
				const auto     to_phi_a = geom::convert_angle_type<angle_base_type>( to_a_res.get_phi_angle() );
				const auto     to_phi_b = geom::convert_angle_type<angle_base_type>( to_b_res.get_phi_angle() );
				return wrapped_difference( to_phi_a, to_phi_b );
			}

			/// \brief Get the (wrapped) difference between the two pairs' to_residue psi angles
			inline angle_type to_psi_angle_difference(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                                          const size_size_pair &prm_indices_b, ///< The indices of the from/to residues in the second protein
			                                          const protein        &prm_protein_a, ///< The first  protein
			                                          const protein        &prm_protein_b  ///< The second protein
			                                          ) {
				const residue    &to_a_res = prm_protein_a.get_residue_ref_of_index( prm_indices_a.second );
				const residue    &to_b_res = prm_protein_b.get_residue_ref_of_index( prm_indices_b.second );
				const angle_type &to_psi_a = geom::make_angle_from_radians<angle_base_type>( angle_in_radians( to_a_res.get_psi_angle() ) );
				const angle_type &to_psi_b = geom::make_angle_from_radians<angle_base_type>( angle_in_radians( to_b_res.get_psi_angle() ) );
				return wrapped_difference( to_psi_a, to_psi_b );
			}

			/// \brief Get the maximum (wrapped) difference between the two pairs' from/to residue phi angles
			inline angle_type max_phi_angle_difference(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                                           const size_size_pair &prm_indices_b, ///< The indices of the from/to residues in the second protein
			                                           const protein        &prm_protein_a, ///< The first  protein
			                                           const protein        &prm_protein_b  ///< The second protein
			                                           ) {
				return std::max(
					from_phi_angle_difference( prm_indices_a, prm_indices_b, prm_protein_a, prm_protein_b ),
					  to_phi_angle_difference( prm_indices_a, prm_indices_b, prm_protein_a, prm_protein_b )
				);
			}

			/// \brief Get the maximum (wrapped) difference between the two pairs' from/to residue psi angles
			inline angle_type max_psi_angle_difference(const size_size_pair &prm_indices_a, ///< The indices of the from/to residues in the first  protein
			                                           const size_size_pair &prm_indices_b, ///< The indices of the from/to residues in the second protein
			                                           const protein        &prm_protein_a, ///< The first  protein
			                                           const protein        &prm_protein_b  ///< The second protein
			                                           ) {
				return std::max(
					from_psi_angle_difference( prm_indices_a, prm_indices_b, prm_protein_a, prm_protein_b ),
					  to_psi_angle_difference( prm_indices_a, prm_indices_b, prm_protein_a, prm_protein_b )
				);
			}

		} // namespace detail
	} // namespace index
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_VCIE_PROTEIN_FUNCTIONS_HPP
