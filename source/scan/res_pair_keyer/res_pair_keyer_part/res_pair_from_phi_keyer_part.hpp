/// \file
/// \brief The res_pair_from_phi_keyer_part class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_FROM_PHI_KEYER_PART_H
#define _CATH_TOOLS_SOURCE_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_FROM_PHI_KEYER_PART_H

//#include <boost/range/irange.hpp>
//#include <boost/range/join.hpp>

#include "scan/detail/res_pair/multi_struc_res_rep_pair.hpp"
#include "scan/quad_criteria.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/detail/res_pair_phi_psi_angle_keyer_part.hpp"
//#include "structure/geometry/angle.hpp"

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Specification for a res_pair_phi_psi_angle_keyer_part for the from_phi angle
			class res_pair_from_phi_keyer_part_spec final {
			public:
				/// \brief Sanity check the specified cell width
				static angle_type sanity_check_cell_width(const angle_type &arg_cell_width ///< The cell width to be sanity-checked
				                                          ) {
					/// \todo Create a `bool is_shifted(const angle &, ...)` helper function for angle and use it here
					return ( arg_cell_width <= geom::zero_angle<angle_base_type>() || arg_cell_width > geom::one_revolution<angle_base_type>() )
						? throw std::logic_error( "Cannot create an angle-based res_pair keyer_part with a cell_width that isn't in ( 0, 2pi ]" )
						: arg_cell_width;
				}

				/// \brief Get a short name that describes this key part
				static std::string get_name() {
					return "from_phi";
				}

				/// \brief Extract the relevant value from the specified res_pair
				static angle_type get_value(const multi_struc_res_rep_pair &arg_res_pair ///< The res_pair from which the relevant value should be extracted
				                            ) {
					return arg_res_pair.get_res_pair_core().get_from_phi_angle();
				}

				/// \brief Extract the search radius from the specified quad_criteria
				static angle_type get_search_radius(const quad_criteria &arg_criteria   ///< The criteria defining what is considered a match
				                                    ) {
					return arg_criteria.get_maximum_phi_angle_difference();
				}
			};

		} // namespace detail

		/// \brief Type alias for keyer_part for from-phi angle
		using res_pair_from_phi_keyer_part = detail::res_pair_phi_psi_angle_keyer_part<detail::res_pair_from_phi_keyer_part_spec>;
	} // namespace scan
} // namespace cath

#endif
