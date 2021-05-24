/// \file
/// \brief The dssp_accessibility class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_ACCESSIBILITY_CALC_DSSP_ACCESSIBILITY_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_ACCESSIBILITY_CALC_DSSP_ACCESSIBILITY_HPP

#include "cath/structure/structure_type_aliases.hpp"

#include <utility>

// clang-format off
namespace cath::file { class pdb; }
namespace cath::file { class pdb_atom; }
namespace cath::file { class pdb_residue; }
namespace cath::file { enum class coarse_element_type : char; }
// clang-format on

namespace cath::sec {

	namespace detail {

		/// \brief Hold a bunch of constants for use in DSSP accessibility calculations
		struct dssp_ball_constants final {

			/// \brief The defalt number with which to specify the sphere of points for accessibility calculations
			static constexpr size_t NUMBER           = 200;

			/// \brief The radius to use in DSSP accessibility calculations for a N     atom
			static constexpr double RADIUS_N         = 1.65;

			/// \brief The radius to use in DSSP accessibility calculations for a CA    atom
			static constexpr double RADIUS_CA        = 1.87;

			/// \brief The radius to use in DSSP accessibility calculations for a C     atom
			static constexpr double RADIUS_C         = 1.76;

			/// \brief The radius to use in DSSP accessibility calculations for a O     atom
			static constexpr double RADIUS_O         = 1.40;

			/// \brief The radius to use in DSSP accessibility calculations for a side  atom
			static constexpr double RADIUS_SIDE_ATOM = 1.80;

			/// \brief The radius to use in DSSP accessibility calculations for a water atom
			static constexpr double RADIUS_WATER     = 1.40;

			/// \brief The maximum distance between two atoms that could possible affect each other's accessibilities
			static constexpr double MAX_ATOM_DIST    = RADIUS_CA + RADIUS_CA + RADIUS_WATER + RADIUS_WATER;
		};

	} // namespace detail

	geom::coord_vec make_dssp_ball_points(const size_t & = detail::dssp_ball_constants::NUMBER);

	double get_dssp_access_radius_without_water(const file::coarse_element_type &);

	double get_dssp_access_radius_without_water(const file::pdb_atom &);

	double get_dssp_access_radius_with_water(const file::coarse_element_type &);

	double get_dssp_access_radius_with_water(const file::pdb_atom &);

	bool access_overlap(const file::pdb_atom &,
	                    const geom::coord &,
	                    const file::pdb_atom &);

	size_t get_accessibility_count(const file::pdb_atom &,
	                               const file::pdb &,
	                               const size_t & = detail::dssp_ball_constants::NUMBER);

	double get_accessibility_fraction(const file::pdb_atom &,
	                                  const file::pdb &,
	                                  const size_t & = detail::dssp_ball_constants::NUMBER);

	double get_accessibility_surface_area(const file::pdb_atom &,
	                                      const file::pdb &,
	                                      const size_t & = detail::dssp_ball_constants::NUMBER);

	double get_accessibility_surface_area(const file::pdb_residue &,
	                                      const file::pdb &,
	                                      const size_t & = detail::dssp_ball_constants::NUMBER);

	doub_vec calc_accessibilities(const file::pdb &,
	                              const size_t & = detail::dssp_ball_constants::NUMBER);

	doub_vec calc_accessibilities_with_scanning(const file::pdb &);

} // namespace cath::sec

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_ACCESSIBILITY_CALC_DSSP_ACCESSIBILITY_HPP
