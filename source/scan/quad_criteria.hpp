/// \file
/// \brief The quad_criteria class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_QUAD_CRITERIA_H
#define _CATH_TOOLS_SOURCE_SCAN_QUAD_CRITERIA_H

#include "scan/detail/scan_type_aliases.hpp"
#include "scan/res_pair_index_dirn_criterion.hpp"
#include "structure/geometry/angle.hpp"

#include <iosfwd>

namespace cath {
	namespace scan {

			/// \brief Define whether a quad (ie a pair of res_pairs) should be counted as a match
			///
			/// At present, this is kept simpler by being read-only after construction.
			///
			/// The criteria are not used for filtering out individual res_pairs in building the scan data structures (only self-self res_pairs are excluded).
			/// That said, the criteria are used to determine the multiple close keys under which each res_pair should be stored in the dense index.
			///  * the final responsibility for filtering matches on the criteria remains with the scanning code
			///    (ie checking every res_pair meets are_not_violated_by() AND every quad meets are_met_by())
			///  * This means that the index must be rebuilt before different criteria can be tried.
			///
			/// All of the individual criteria must all be satisfied for the quad_criteria to be satisfied
			class quad_criteria final {
			private:
				/// \brief Whether the two res_pairs are required to have matching directions
				///       (ie either both have from_index < to_index or both have from_index > to_index)
				res_pair_index_dirn_criterion index_direction_criterion;

				/// \brief The minimum distance required between each res_pair's from_index and its to_index
				///
				/// This criterion is applied to a single res_pair, rather than a pair, so a function operator
				/// is provided to perform that test
				index_type minimum_index_distance;

				/// \brief The maximum squared distance permissible between the two res_pairs' views
				detail::view_base_type maximum_squared_distance;

				/// \brief The maximum distance 1 permissible between the two res_pairs' view frame quaternions
				///
				/// This is exactly equivalent to the maximum_frame_angle_difference test below and is populated automatically
				/// based on the angle provided for that but it can be calculated a bit quicker. Quaternions with no angle between
				/// them give a distance of 0 and quaternions at 180 degrees from each other give a distance of 1.
				detail::frame_quat_rot_type maximum_frame_angle_distance_1;

				/// \brief The maximum angle permissible between the two res_pairs' view frames
				detail::angle_type maximum_frame_angle_difference;

				/// \brief The maximum angle permissible between the two res_pairs' phi angles
				///
				/// This should handle wrapping so that similar angles should always be recognised as such
				detail::angle_type maximum_phi_angle_difference;

				/// \brief The maximum angle permissible between the two res_pairs' psi angles
				///
				/// This should handle wrapping so that similar angles should always be recognised as such
				detail::angle_type maximum_psi_angle_difference;

			public:
				quad_criteria(const res_pair_index_dirn_criterion &,
				              const index_type &,
				              const detail::view_base_type &,
				              const detail::angle_type &,
				              const detail::angle_type &,
				              const detail::angle_type &);

				const res_pair_index_dirn_criterion & get_index_direction_criterion() const;
				const index_type & get_minimum_index_distance() const;
				const detail::view_base_type & get_maximum_squared_distance() const;
				const detail::frame_quat_rot_type & get_maximum_frame_angle_distance_1() const;
				const detail::angle_type & get_maximum_frame_angle_difference() const;
				const detail::angle_type & get_maximum_phi_angle_difference() const;
				const detail::angle_type & get_maximum_psi_angle_difference() const;
			};

			/// \brief Getter for require_matching_directions
			inline const res_pair_index_dirn_criterion & quad_criteria::get_index_direction_criterion() const {
				return index_direction_criterion;
			}

			/// \brief Getter for minimum_index_distance
			inline const index_type & quad_criteria::get_minimum_index_distance() const {
				return minimum_index_distance;
			}

			/// \brief Getter for maximum_squared_distance
			inline const detail::view_base_type & quad_criteria::get_maximum_squared_distance() const {
				return maximum_squared_distance;
			}

			/// \brief Getter for maximum_frame_angle_distance_1
			inline const detail::frame_quat_rot_type & quad_criteria::get_maximum_frame_angle_distance_1() const {
				return maximum_frame_angle_distance_1;
			}

			/// \brief Getter for maximum_frame_angle_difference
			inline const detail::angle_type & quad_criteria::get_maximum_frame_angle_difference() const {
				return maximum_frame_angle_difference;
			}

			/// \brief Getter for maximum_phi_angle_difference
			inline const detail::angle_type & quad_criteria::get_maximum_phi_angle_difference() const {
				return maximum_phi_angle_difference;
			}

			/// \brief Getter for maximum_psi_angle_difference
			inline const detail::angle_type & quad_criteria::get_maximum_psi_angle_difference() const {
				return maximum_psi_angle_difference;
			}

			/// \brief Getter for require_matching_directions
			inline bool requires_matching_directions(const quad_criteria &arg_criteria ///< TODOCUMENT
			                                         ) {
				return arg_criteria.get_index_direction_criterion() == res_pair_index_dirn_criterion::MUST_MATCH;
			}

			/// \brief Get the maximum allowed distance implied by the specified quad_criteria
			inline detail::view_base_type get_maximum_distance(const quad_criteria &arg_quad_criteria ///< The arg_quad_criteria to query
			                                                   ) {
				return std::sqrt( arg_quad_criteria.get_maximum_squared_distance() );
			}

			quad_criteria make_default_quad_criteria();
			quad_criteria parse_quad_criteria(const std::string &);
			quad_criteria_vec get_standard_quad_criterias();

			std::ostream & operator<<(std::ostream &,
			                          const quad_criteria &);
	} // namespace scan
} // namespace cath

#endif
