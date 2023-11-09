/// \file
/// \brief The view_cache_index_entry class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_VIEW_CACHE_INDEX_ENTRY_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_VIEW_CACHE_INDEX_ENTRY_HPP

#include <iostream> // To fix "boost/geometry/algorithms/detail/overlay/handle_colocations.hpp:198:10: error: no member named 'cout' in namespace 'std'" on MacOS travis-ci build with Boost 1.61.0

#include <boost/geometry/algorithms/comparable_distance.hpp>

#include "cath/common/difference.hpp" /// ***** TEMPORARY *****
#include "cath/common/exception/invalid_argument_exception.hpp" /// ***** TEMPORARY *****
#include "cath/ssap/context_res.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/quat_rot.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/view_cache/index/detail/vcie_protein_functions.hpp"
#include "cath/structure/view_cache/index/detail/view_cache_index_type_aliases.hpp"

#include <utility>

namespace cath::index {

	/// \brief Entry describing a pair of from/to residues in way that allows for quick searching for
	///        similar pairs of such entries in a view_cache_index
	///
	/// Invariants:
	///  * All four from/to phi/psi angles must be shifted to the same (standard) angle range
	class view_cache_index_entry final {
	private:
		/// \brief The view of the to_residue from the from_residue
		/// (ie the vector from the from_residue to the to_residue, in the coordinate frame of the from_residue)
		detail::view_type view;

		/// \brief The coordinate frame of the from_residue,
		///        as determined by its core atoms
		detail::frame_quat_rot frame;

		/// \brief The phi angle of the from_residue
		detail::angle_type from_phi_angle;

		/// \brief The psi angle of the from_residue
		detail::angle_type from_psi_angle;

		/// \brief The phi angle of the to_residue
		detail::angle_type to_phi_angle;

		/// \brief The psi angle of the to_residue
		detail::angle_type to_psi_angle;

					/// \brief The index of the from_residue in its protein
		detail::index_type from_index;

		/// \brief The index of the to_residue in its protein
		detail::index_type to_index;

	public:
		view_cache_index_entry(const detail::index_type &,
		                       const detail::index_type &,
		                       detail::view_type,
		                       const geom::rotation &,
		                       const detail::angle_type &,
		                       const detail::angle_type &,
		                       const detail::angle_type &,
		                       const detail::angle_type &);

		[[nodiscard]] const detail::index_type &    get_from_index() const;
		[[nodiscard]] const detail::index_type &    get_to_index() const;
		[[nodiscard]] const detail::view_type &     get_view() const;
		[[nodiscard]] const detail::frame_quat_rot &get_frame() const;
		[[nodiscard]] const detail::angle_type &    get_from_phi_angle() const;
		[[nodiscard]] const detail::angle_type &    get_from_psi_angle() const;
		[[nodiscard]] const detail::angle_type &    get_to_phi_angle() const;
		[[nodiscard]] const detail::angle_type &    get_to_psi_angle() const;
	};

	/// \brief Getter for from_index
	inline const detail::index_type & view_cache_index_entry::get_from_index() const {
		return from_index;
	}

	/// \brief Getter for to_index
	inline const detail::index_type & view_cache_index_entry::get_to_index() const {
		return to_index;
	}

	/// \brief Getter for view
	inline const detail::view_type & view_cache_index_entry::get_view() const {
		return view;
	}

	/// \brief Getter for frame
	inline const detail::frame_quat_rot & view_cache_index_entry::get_frame() const {
		return frame;
	}

	/// \brief Getter for from_phi_angle
	inline const detail::angle_type & view_cache_index_entry::get_from_phi_angle() const {
		return from_phi_angle;
	}

	/// \brief Getter for from_psi_angle
	inline const detail::angle_type & view_cache_index_entry::get_from_psi_angle() const {
		return from_psi_angle;
	}

	/// \brief Getter for to_phi_angle
	inline const detail::angle_type & view_cache_index_entry::get_to_phi_angle() const {
		return to_phi_angle;
	}

	/// \brief Getter for to_psi_angle
	inline const detail::angle_type & view_cache_index_entry::get_to_psi_angle() const {
		return to_psi_angle;
	}

	std::ostream & operator<<(std::ostream &,
	                          const view_cache_index_entry &);

	namespace detail {

		view_cache_index_entry make_view_cache_index_entry(const protein &,
		                                                   const size_t &,
		                                                   const size_t &);


		/// \brief Return the absolute difference between the from_index and to_index of the specified view_cache_index_entry
		inline detail::index_type index_difference(const view_cache_index_entry &prm_cache ///< The view_cache_index_entry to query
		                                           ) {
			return common::difference( prm_cache.get_from_index(), prm_cache.get_to_index() );
		}

		/// \brief Calculate the distance between the quaternions for the two pairs' view_frames
		///        (the coordinate frames of their to_residues as seen from the coordinate frames of their from_residues)
		///
		/// This distance can be calculated a bit quicker than the angle and can be used in an exactly equivalent criterion
		///
		/// \relates view_cache_index_entry
		inline detail::frame_quat_rot_type distance_1_between_frames(const view_cache_index_entry &prm_cache_a, ///< The first  view_cache_index_entry
		                                                             const view_cache_index_entry &prm_cache_b  ///< The second view_cache_index_entry
		                                                             ) {
			return distance_1_between_quat_rots( prm_cache_a.get_frame(), prm_cache_b.get_frame() );
		}

		/// \brief Calculate the squared distance between the views of the two residue pairs
		///
		/// Each view is the location of the to_residue's carbon-beta atom as seen from
		/// the coordinate frame of the from_residue
		///
		/// \relates view_cache_index_entry
		inline double squared_distance(const view_cache_index_entry &prm_cache_a, ///< The first  view_cache_index_entry
		                               const view_cache_index_entry &prm_cache_b  ///< The second view_cache_index_entry
		                               ) {
			// const auto the_distance = boost::geometry::distance(
			// 	prm_cache_a.get_view(),
			// 	prm_cache_b.get_view()
			// );
			// const auto the_comp_distance = boost::geometry::comparable_distance(
			// 	prm_cache_a.get_view(),
			// 	prm_cache_b.get_view()
			// );
			// if ( common::difference( the_comp_distance, the_distance * the_distance ) > 0.0001 ) {
			// 	std::cerr << "the_distance is "                << the_distance                << std::endl;
			// 	std::cerr << "the_distance * the_distance is " << the_distance * the_distance << std::endl;
			// 	std::cerr << "the_comp_distance is "           << the_comp_distance           << std::endl;
			// 	BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Distances do not match"));
			// }
			return boost::geometry::comparable_distance(
				prm_cache_a.get_view(),
				prm_cache_b.get_view()
			);
			// return squared_distance_between_points(
			// 	prm_cache_a.get_view(),
			// 	prm_cache_b.get_view()
			// );
		}

		/// \brief Return whether the two pairs of residues are both in the same direction
		///        (ie both have from_index < to_index or both have from_index > to_index)
		///
		/// \relates view_cache_index_entry
		inline bool same_direction(const view_cache_index_entry &prm_cache_a, ///< The first  view_cache_index_entry
		                           const view_cache_index_entry &prm_cache_b  ///< The second view_cache_index_entry
		                           ) {
			return (
				( prm_cache_a.get_from_index() < prm_cache_a.get_to_index() )
				==
				( prm_cache_b.get_from_index() < prm_cache_b.get_to_index() )
			);
		}

		/// \brief Calculate the minimum of the two pairs' absolute differences between their from_index and their to_index
		///
		/// \relates view_cache_index_entry
		inline size_t min_index_difference(const view_cache_index_entry &prm_cache_a, ///< The first  view_cache_index_entry
		                                   const view_cache_index_entry &prm_cache_b  ///< The second view_cache_index_entry
		                                   ) {
			return min_index_difference(
				std::make_pair( prm_cache_a.get_from_index(), prm_cache_a.get_to_index() ),
				std::make_pair( prm_cache_b.get_from_index(), prm_cache_b.get_to_index() )
			);
		}

		/// \brief Get the (wrapped) difference between the two pairs' from_residue phi angles
		///
		/// \relates view_cache_index_entry
		inline detail::angle_type from_phi_angle_difference(const view_cache_index_entry &prm_cache_a, ///< The first  view_cache_index_entry
		                                                    const view_cache_index_entry &prm_cache_b  ///< The second view_cache_index_entry
		                                                    ) {
			return unshifted_wrapped_difference(
				prm_cache_a.get_from_phi_angle(),
				prm_cache_b.get_from_phi_angle()
			);
		}

		/// \brief Get the (wrapped) difference between the two pairs' from_residue psi angles
		///
		/// \relates view_cache_index_entry
		inline detail::angle_type from_psi_angle_difference(const view_cache_index_entry &prm_cache_a, ///< The first  view_cache_index_entry
		                                                    const view_cache_index_entry &prm_cache_b  ///< The second view_cache_index_entry
		                                                    ) {
			return unshifted_wrapped_difference(
				prm_cache_a.get_from_psi_angle(),
				prm_cache_b.get_from_psi_angle()
			);
		}

		/// \brief Get the (wrapped) difference between the two pairs' to_residue phi angles
		///
		/// \relates view_cache_index_entry
		inline detail::angle_type to_phi_angle_difference(const view_cache_index_entry &prm_cache_a,   ///< The first  view_cache_index_entry
		                                                  const view_cache_index_entry &prm_cache_b    ///< The second view_cache_index_entry
		                                                  ) {
			return unshifted_wrapped_difference(
				prm_cache_a.get_to_phi_angle(),
				prm_cache_b.get_to_phi_angle()
			);
		}

		/// \brief Get the (wrapped) difference between the two pairs' to_residue psi angles
		///
		/// \relates view_cache_index_entry
		inline detail::angle_type to_psi_angle_difference(const view_cache_index_entry &prm_cache_a,   ///< The first  view_cache_index_entry
		                                                  const view_cache_index_entry &prm_cache_b    ///< The second view_cache_index_entry
		                                                  ) {
			return unshifted_wrapped_difference(
				prm_cache_a.get_to_psi_angle(),
				prm_cache_b.get_to_psi_angle()
			);
		}

		/// \brief Get the maximum (wrapped) difference between the two pairs' from/to residue phi angles
		///
		/// \relates view_cache_index_entry
		inline detail::angle_type max_phi_angle_difference(const view_cache_index_entry &prm_cache_a,   ///< The first  view_cache_index_entry
		                                                   const view_cache_index_entry &prm_cache_b    ///< The second view_cache_index_entry
		                                                   ) {
			return std::max(
				from_phi_angle_difference( prm_cache_a, prm_cache_b ),
				  to_phi_angle_difference( prm_cache_a, prm_cache_b )
			);
		}

		/// \brief Get the maximum (wrapped) difference between the two pairs' from/to residue psi angles
		///
		/// \relates view_cache_index_entry
		inline detail::angle_type max_psi_angle_difference(const view_cache_index_entry &prm_cache_a,   ///< The first  view_cache_index_entry
		                                                   const view_cache_index_entry &prm_cache_b    ///< The second view_cache_index_entry
		                                                   ) {
			return std::max(
				from_psi_angle_difference( prm_cache_a, prm_cache_b ),
				  to_psi_angle_difference( prm_cache_a, prm_cache_b )
			);
		}

		/// \brief Convenience function to get the x component of the view in the specified view_cache_index_entry
		///
		/// \relates view_cache_index_entry
		inline const view_base_type & get_view_x(const view_cache_index_entry &prm_cache ///< The view_cache_index_entry to query
		                                          ) {
			return prm_cache.get_view().get<0>();
		}

		/// \brief Convenience function to get the y component of the view in the specified view_cache_index_entry
		///
		/// \relates view_cache_index_entry
		inline const view_base_type & get_view_y(const view_cache_index_entry &prm_cache ///< The view_cache_index_entry to query
		                                          ) {
			return prm_cache.get_view().get<1>();
		}

		/// \brief Convenience function to get the z component of the view in the specified view_cache_index_entry
		///
		/// \relates view_cache_index_entry
		inline const view_base_type & get_view_z(const view_cache_index_entry &prm_cache ///< The view_cache_index_entry to query
		                                          ) {
			return prm_cache.get_view().get<2>();
		}

	} // namespace detail

} // namespace cath::index

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_VIEW_CACHE_INDEX_ENTRY_HPP
