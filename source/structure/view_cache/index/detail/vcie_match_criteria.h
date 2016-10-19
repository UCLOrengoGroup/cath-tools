/// \file
/// \brief The vcie_match_criteria class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_VCIE_MATCH_CRITERIA_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_VCIE_MATCH_CRITERIA_H

#include "structure/geometry/angle.h"
#include "structure/view_cache/index/view_cache_index_entry.h"

namespace cath { class protein; }
namespace cath { namespace index { class view_cache_index_entry; } }

/// \todo Consider moving this outside of the detail sub-namespace because users might
///       well need to interact with vcie_match_criteria. In particular, it doesn't
///       make much sense to hide make_default_vcie_match_criteria() away in detail::

namespace cath {
	namespace index {
		namespace detail {

			/// \brief Store the criteria for determining whether a pair of view_cache_index_entries (vcies)
			///        merit further consideration
			///
			/// A view_cache_index_entry captures information about a pair of residues in a protein (one "from" and one "to)
			/// and this criteria class determines whether it's worth considering one from-to pair in one protein structure with
			/// another from-to pair in another protein structure.
			///
			/// At present, this is kept simpler by being read-only after construction.
			///
			/// All of the individual criteria must all be satisfied for the vcie_match_criteria to be satisfied
			class vcie_match_criteria final {
			private:
				/// \brief Whether the two vcies are required to have matching directions
				///       (ie either both have from_index < to_index or both have from_index > to_index)
				bool require_matching_directions;

				/// \brief The minimum distance required between each vcie's from_index and its to_index
				///
				/// This criterion is applied to a single vcie, rather than a pair, so a function operator
				/// is provided to perform that test
				index_type minimum_index_distance;

				/// \brief The maximum squared distance permissible between the two vcies' views
				view_base_type maximum_squared_distance;

				/// \brief The maximum distance 1 permissible between the two vcies' view frame quaternions
				///
				/// This is exactly equivalent to the maximum_frame_angle_difference test below and is populated automatically
				/// based on the angle provided for that but it can be calculated a bit quicker. Quaternions with no angle between
				/// them give a distance of 0 and quaternions at 180 degrees from each other give a distance of 1.
				frame_quat_rot_type maximum_frame_angle_distance_1;

				/// \brief The maximum angle permissible between the two vcies' view frames
				angle_type maximum_frame_angle_difference;

				/// \brief The maximum angle permissible between the two vcies' phi angles
				///
				/// This should handle wrapping so that similar angles should always be recognised as such
				angle_type maximum_phi_angle_difference;

				/// \brief The maximum angle permissible between the two vcies' psi angles
				///
				/// This should handle wrapping so that similar angles should always be recognised as such
				angle_type maximum_psi_angle_difference;

			public:
				vcie_match_criteria(const bool &,
				                    const index_type &,
				                    const view_base_type &,
				                    const angle_type &,
				                    const angle_type &,
				                    const angle_type &);

				const bool & get_require_matching_directions() const;
				const index_type & get_minimum_index_distance() const;
				const view_base_type & get_maximum_squared_distance() const;
				const frame_quat_rot_type & get_maximum_frame_angle_distance_1() const;
				const angle_type & get_maximum_frame_angle_difference() const;
				const angle_type & get_maximum_phi_angle_difference() const;
				const angle_type & get_maximum_psi_angle_difference() const;

				bool operator()(const view_cache_index_entry &) const;

//				template <bool CHECK_DIRN = false>
				bool operator()(const view_cache_index_entry &,
				                const view_cache_index_entry &) const;

				bool operator()(const size_size_pair &,
				                const size_size_pair &,
				                const protein &,
				                const protein &) const;
			};

			/// \brief Getter for require_matching_directions
			inline const bool & vcie_match_criteria::get_require_matching_directions() const {
				return require_matching_directions;
			}

			/// \brief Getter for minimum_index_distance
			inline const index_type & vcie_match_criteria::get_minimum_index_distance() const {
				return minimum_index_distance;
			}

			/// \brief Getter for maximum_squared_distance
			inline const view_base_type & vcie_match_criteria::get_maximum_squared_distance() const {
				return maximum_squared_distance;
			}

			/// \brief Getter for maximum_frame_angle_distance_1
			inline const frame_quat_rot_type & vcie_match_criteria::get_maximum_frame_angle_distance_1() const {
				return maximum_frame_angle_distance_1;
			}

			/// \brief Getter for maximum_frame_angle_difference
			inline const angle_type & vcie_match_criteria::get_maximum_frame_angle_difference() const {
				return maximum_frame_angle_difference;
			}

			/// \brief Getter for maximum_phi_angle_difference
			inline const angle_type & vcie_match_criteria::get_maximum_phi_angle_difference() const {
				return maximum_phi_angle_difference;
			}

			/// \brief Getter for maximum_psi_angle_difference
			inline const angle_type & vcie_match_criteria::get_maximum_psi_angle_difference() const {
				return maximum_psi_angle_difference;
	}

			/// \brief Returns whether an individual view_cache_index_entry has any chance of matching another
			///
			/// This method allows some view_cache_index_entries to be ruled out of consideration for matching
			/// without having it to compare it to others.
			///
			/// This can be the case if the from_index matches the to_index or their absolute difference is
			/// otherwise less than the criteria's minimum_index_distance
			inline bool vcie_match_criteria::operator()(const view_cache_index_entry &arg_cache ///< The view_cache_index_entry to be tested for whether it has any chance of matching another
			                                            ) const {
				if ( arg_cache.get_from_index()    == arg_cache.get_to_index()     ) {
					return false;
				}
				if ( index_difference( arg_cache ) <  get_minimum_index_distance() ) {
					return false;
				}
				return true;
			}

			/// \brief Return whether the two specified vcies match (if they both pass the single vcie test, above)
			///
			/// \todo The require_matching_directions test has been removed because this function needs to be very fast and
			///       the view_cache_index should already ensure that only matching direction entries are compared.
			///       However, it may cause confusion that the function can can now accept naked pairs of entries with different directions.
			///       Try to reduce the chance of that unexpected behaviour biting some poor unsuspecting soul.
			///
			/// IMPORTANT: This result should only be considered if both vcies also pass the single vcie test above.
//			template <bool CHECK_DIRN = false>
			inline bool vcie_match_criteria::operator()(const view_cache_index_entry &arg_cache_a,   ///< The first  view_cache_index_entry to test
			                                            const view_cache_index_entry &arg_cache_b    ///< The second view_cache_index_entry to test
			                                            ) const {
//				if ( CHECK_DIRN ) {
//					if ( get_require_matching_directions() && ! same_direction( arg_indices_a, arg_indices_b ) ) {
//						return false;
//					}
//				}
				if ( squared_distance         ( arg_cache_a, arg_cache_b ) > get_maximum_squared_distance()       ) {
					return false;
				}
				if ( distance_1_between_frames( arg_cache_a, arg_cache_b ) > get_maximum_frame_angle_distance_1() ) {
					return false;
				}
				if ( max_phi_angle_difference ( arg_cache_a, arg_cache_b ) > get_maximum_phi_angle_difference()   ) {
					return false;
				}
				if ( max_psi_angle_difference ( arg_cache_a, arg_cache_b ) > get_maximum_psi_angle_difference()   ) {
					return false;
				}
				return true;
			}

			/// \brief Return whether a quad of residues are considered a match
			///        by using only the raw indices and proteins
			inline bool vcie_match_criteria::operator()(const size_size_pair &arg_indices_a, ///< The from_index and from_index for the first  protein
			                                            const size_size_pair &arg_indices_b, ///< The from_index and from_index for the second protein
			                                            const protein        &arg_protein_a, ///< The first  protein
			                                            const protein        &arg_protein_b  ///< The second protein
			                                            ) const {
// #ifndef NDEBUG
// 				// Until the class's internal protein references get removed, check that they refer to the
// 				// same objects as the ones passed as arguments here
// 				if ( &arg_protein_a != &( protein_a.get() ) ) {
// 					BOOST_THROW_EXCEPTION(cath::common::out_of_range_exception("vcie_match_criteria function operator called with different protein_a from the one used to initialise the vcie_match_criteria object"));
// 				}
// 				if ( &arg_protein_b != &( protein_b.get() ) ) {
// 					BOOST_THROW_EXCEPTION(cath::common::out_of_range_exception("vcie_match_criteria function operator called with different protein_b from the one used to initialise the vcie_match_criteria object"));
// 				}
// #endif
				if ( get_require_matching_directions() && ! same_direction( arg_indices_a, arg_indices_b ) ) {
					return false;
				}
				if ( min_index_difference      ( arg_indices_a, arg_indices_b                               ) < get_minimum_index_distance()         ) {
					return false;
				}
				if ( squared_distance          ( arg_indices_a, arg_indices_b, arg_protein_a, arg_protein_b ) > get_maximum_squared_distance()       ) {
					return false;
				}
				if ( angle_between_frames      ( arg_indices_a, arg_indices_b, arg_protein_a, arg_protein_b ) > get_maximum_frame_angle_difference() ) {
					return false;
				}
				if ( max_phi_angle_difference  ( arg_indices_a, arg_indices_b, arg_protein_a, arg_protein_b ) > get_maximum_phi_angle_difference()   ) {
					return false;
				}
				if ( max_psi_angle_difference  ( arg_indices_a, arg_indices_b, arg_protein_a, arg_protein_b ) > get_maximum_psi_angle_difference()   ) {
					return false;
				}
				return true;
			}

			vcie_match_criteria make_default_vcie_match_criteria();
			vcie_match_criteria parse_vcie_match_criteria(const std::string &);
			vcie_match_criteria_vec get_standard_vcie_match_criterias();

			std::ostream & operator<<(std::ostream &,
			                          const vcie_match_criteria &);
		} // namespace detail
	} // namespace index
} // namespace cath

#endif
