/// \file
/// \brief The res_pair_view_y_keyer_part class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_VIEW_Y_KEYER_PART_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_VIEW_Y_KEYER_PART_HPP

#include "cath/scan/detail/res_pair/multi_struc_res_rep_pair.hpp"
#include "cath/scan/detail/res_pair/res_pair_core.hpp"
#include "cath/scan/quad_criteria.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer_part/detail/axis_keyer_part.hpp"

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Specification for a res_pair_view_axis_keyer_part for the view's y-dimension
			template <typename Stored, typename Crit>
			class res_pair_view_y_keyer_part_spec final {
			public:
				/// \brief TODOCUMENT
				using stored_t   = Stored;

				/// \brief TODOCUMENT
				using criteria_t = Crit;

				/// \brief TODOCUMENT
				using value_t    = view_base_type;

				/// \brief Sanity check the specified cell width
				static constexpr value_t sanity_check_cell_width(const value_t &prm_cell_width ///< The cell width to be sanity-checked
				                                                 ) {
					return ( prm_cell_width <= 0.0 || prm_cell_width > 16384.0 ) ? throw std::logic_error( "Cannot create an axis-based res_pair keyer_part with a cell_width that isn't a sensible, strictly-positive value" )
					                                                             : prm_cell_width;
				}

				/// \brief Get a short name that describes this key part
				static std::string get_name() {
					return "view_y";
				}

				/// \brief Extract the relevant value from the specified res_pair
				static constexpr value_t get_value(const stored_t &prm_res_pair ///< The res_pair from which the relevant value should be extracted
				                                   ) {
					return get_view_y( prm_res_pair );
				}

				/// \brief Extract the search radius from the specified quad_criteria
				static constexpr value_t get_search_radius(const criteria_t &prm_criteria ///< The criteria defining what is considered a match
				                                           ) {
					return get_maximum_distance( prm_criteria );
				}
			};

		} // namespace detail

		/// \brief Type alias for keyer_part for y-axis of view
		using res_pair_view_y_keyer_part = detail::axis_keyer_part<detail::res_pair_view_y_keyer_part_spec< detail::multi_struc_res_rep_pair, quad_criteria > >;
	} // namespace scan
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_VIEW_Y_KEYER_PART_HPP
