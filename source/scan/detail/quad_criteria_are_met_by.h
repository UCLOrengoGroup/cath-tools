/// \file
/// \brief The quad_criteria are_met_by functions header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_DETAIL_QUAD_CRITERIA_ARE_MET_BY_H
#define _CATH_TOOLS_SOURCE_SCAN_DETAIL_QUAD_CRITERIA_ARE_MET_BY_H

#include "scan/detail/res_pair/functions/res_index_pair_functions.h"
#include "scan/detail/res_pair/functions/res_pair_core_functions.h"
#include "scan/detail/res_pair/multi_struc_res_rep_pair.h"
#include "scan/detail/res_pair/single_struc_res_pair.h"
#include "scan/quad_criteria.h"

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Returns whether a the specified criteria are met by the specified single_struc_res_pair
			///
			/// NOTE: As indicated by the name, this can be used to rule some single_struc_res_rep_pairs out of
			///       consideration but the criteria assess the similarity of two res_pairs so it doesn't make sense
			///       to ask whether the criteria are met by a single res_pair.
			///
			/// Violation can occur if the from_index matches the to_index or their absolute difference is
			/// otherwise less than the criteria's minimum_index_distance
			///
			/// \relates quad_criteria
			inline bool are_not_violated_by(const quad_criteria &arg_criteria,       ///< The criteria to apply
			                                const index_type    &arg_from_res_index, ///< The res_pair to be tested for whether it has any chance of matching another
			                                const index_type    &arg_to_res_index    ///< The res_pair to be tested for whether it has any chance of matching another
			                                ) {
				if ( arg_from_res_index == arg_to_res_index ) {
					return false;
				}
				if ( common::difference( arg_from_res_index, arg_to_res_index ) < arg_criteria.get_minimum_index_distance() ) {
					return false;
				}
				return true;
			}

			/// \brief Returns whether a the specified criteria are met by the specified multi_struc_res_pair
			///
			/// NOTE: As indicated by the name, this can be used to rule some single_struc_res_rep_pairs out of
			///       consideration but the criteria assess the similarity of two res_pairs so it doesn't make sense
			///       to ask whether the criteria are met by a single res_pair.
			///
			/// Violation can occur if the from_index matches the to_index or their absolute difference is
			/// otherwise less than the criteria's minimum_index_distance
			///
			/// \relates quad_criteria
			inline bool are_not_violated_by(const quad_criteria            &/*arg_criteria*/, ///< The criteria to apply
			                                const multi_struc_res_rep_pair &arg_res_pair      ///< The res_pair to be tested for whether it has any chance of matching another
			                                ) {
				const auto from_res_rep_index = arg_res_pair.get_from_res_rep_index();
				const auto to_res_rep_index   = arg_res_pair.get_to_res_rep_index  ();
				return ( from_res_rep_index != to_res_rep_index );
			}

			/// \brief Returns whether a the specified criteria are met by the specified single_struc_res_pair
			///
			/// NOTE: As indicated by the name, this can be used to rule some single_struc_res_rep_pairs out of
			///       consideration but the criteria assess the similarity of two res_pairs so it doesn't make sense
			///       to ask whether the criteria are met by a single res_pair.
			///
			/// Violation can occur if the from_index matches the to_index or their absolute difference is
			/// otherwise less than the criteria's minimum_index_distance
			///
			/// \relates quad_criteria
			inline bool are_not_violated_by(const quad_criteria         &arg_criteria, ///< The criteria to apply
			                                const single_struc_res_pair &arg_res_pair  ///< The res_pair to be tested for whether it has any chance of matching another
			                                ) {
				const auto from_res_index = arg_res_pair.get_from_res_idx();
				const auto to_res_index   = arg_res_pair.get_to_res_idx  ();
				if ( from_res_index == to_res_index ) {
					return false;
				}
				if ( common::difference( from_res_index, to_res_index ) < arg_criteria.get_minimum_index_distance() ) {
					return false;
				}
				return true;
			}


			/// \brief Return whether the two specified res_pair_cores match
			///
			/// \relates quad_criteria
			inline bool are_met_by(const quad_criteria &arg_criteria,   ///< The criteria to apply
			                       const res_pair_core &arg_res_pair_a, ///< The first  res_pair to test
			                       const res_pair_core &arg_res_pair_b  ///< The second res_pair to test
			                       ) {
				if ( squared_distance         ( arg_res_pair_a, arg_res_pair_b ) > arg_criteria.get_maximum_squared_distance()       ) {
					return false;
				}
				if ( distance_1_between_frames( arg_res_pair_a, arg_res_pair_b ) > arg_criteria.get_maximum_frame_angle_distance_1() ) {
					return false;
				}
				if ( max_phi_angle_difference ( arg_res_pair_a, arg_res_pair_b ) > arg_criteria.get_maximum_phi_angle_difference()   ) {
					return false;
				}
				if ( max_psi_angle_difference ( arg_res_pair_a, arg_res_pair_b ) > arg_criteria.get_maximum_psi_angle_difference()   ) {
					return false;
				}
				return true;
			}

			/// \brief Return whether the two specified res_pairs match (*assuming they both pass the are_not_violated_by()*)
			///
			/// \todo The require_matching_directions test has been removed because this function needs to be very fast and
			///       the view_cache_index should already ensure that only matching direction entries are compared.
			///       However, it may cause confusion that the function can can now accept naked pairs of entries with different directions.
			///       Try to reduce the chance of that unexpected behaviour biting some poor unsuspecting soul.
			///
			/// IMPORTANT: This result should only be considered if both res_pairs also pass the are_not_violated_by() test above.
			///
			/// \relates quad_criteria
//			template <bool CHECK_DIRN = false>
			inline bool are_met_by(const quad_criteria                    &arg_criteria,   ///< The criteria to apply
			                       const multi_struc_res_rep_pair &arg_res_pair_a, ///< The first  res_pair to test
			                       const multi_struc_res_rep_pair &arg_res_pair_b  ///< The second res_pair to test
			                       ) {
				if ( ! are_met_by( arg_criteria, arg_res_pair_a.get_res_pair_core(), arg_res_pair_b.get_res_pair_core() ) ) {
					return false;
				}
				if ( requires_matching_directions( arg_criteria ) && ! same_direction( arg_res_pair_a, arg_res_pair_b ) ) {
					return false;
				}
				return true;
			}

			/// \brief Return whether the two specified res_pairs match (*assuming they both pass are_not_violated_by()*)
			///
			/// IMPORTANT: This result should only be considered if both res_pairs also pass the are_not_violated_by() test above.
			///
			/// \relates quad_criteria
			inline bool are_met_by(const quad_criteria         &arg_criteria,   ///< The criteria to apply
			                       const single_struc_res_pair &arg_res_pair_a, ///< The first  res_pair to test
			                       const single_struc_res_pair &arg_res_pair_b  ///< The second res_pair to test
			                       ) {
				if ( ! are_met_by( arg_criteria, arg_res_pair_a.get_res_pair_core(), arg_res_pair_b.get_res_pair_core() ) ) {
					return false;
				}
				if ( requires_matching_directions( arg_criteria ) && ! same_direction( arg_res_pair_a, arg_res_pair_b ) ) {
					return false;
				}
				return true;
			}


			/// \brief Return whether a quad of residues is considered a match
			///        by using only the raw indices and proteins
			///
			/// \relates quad_criteria
			inline bool are_met_by(const quad_criteria  &arg_criteria,  ///< The criteria to apply
			                       const size_size_pair &arg_indices_a, ///< The from_index and from_index for the first  protein
			                       const size_size_pair &arg_indices_b, ///< The from_index and from_index for the second protein
			                       const protein        &arg_protein_a, ///< The first  protein
			                       const protein        &arg_protein_b  ///< The second protein
			                       ) {
				if ( requires_matching_directions( arg_criteria ) && ! same_direction( arg_indices_a, arg_indices_b ) ) {
					return false;
				}
				if ( min_index_difference      ( arg_indices_a, arg_indices_b                               ) < arg_criteria.get_minimum_index_distance()         ) {
					return false;
				}
				if ( squared_distance          ( arg_indices_a, arg_indices_b, arg_protein_a, arg_protein_b ) > arg_criteria.get_maximum_squared_distance()       ) {
					return false;
				}
				if ( angle_between_frames      ( arg_indices_a, arg_indices_b, arg_protein_a, arg_protein_b ) > arg_criteria.get_maximum_frame_angle_difference() ) {
					return false;
				}
				if ( max_phi_angle_difference  ( arg_indices_a, arg_indices_b, arg_protein_a, arg_protein_b ) > arg_criteria.get_maximum_phi_angle_difference()   ) {
					return false;
				}
				if ( max_psi_angle_difference  ( arg_indices_a, arg_indices_b, arg_protein_a, arg_protein_b ) > arg_criteria.get_maximum_psi_angle_difference()   ) {
					return false;
				}
				return true;
			}
		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
