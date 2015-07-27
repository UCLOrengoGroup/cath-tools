/// \file
/// \brief The res_pair_view_x_keyer_part class header

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

#ifndef RES_PAIR_VIEW_X_KEYER_PART_H_INCLUDED
#define RES_PAIR_VIEW_X_KEYER_PART_H_INCLUDED

#include "scan/detail/res_pair/multi_struc_res_rep_pair.h"
#include "scan/detail/res_pair/res_pair_core.h"
#include "scan/quad_criteria.h"
#include "scan/res_pair_keyer/res_pair_keyer_part/detail/res_pair_view_axis_keyer_part.h"

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Specification for a res_pair_view_axis_keyer_part for the view's x-dimension
			class res_pair_view_x_keyer_part_spec final {
			public:
				/// \brief Sanity check the specified cell width
				static void sanity_check_cell_width(const view_base_type &arg_cell_width ///< The cell width to be sanity-checked
				                                    ) {
					if ( ! boost::math::isnormal( arg_cell_width ) || arg_cell_width <= 0.0 || arg_cell_width > 16384.0 ) {
						BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot create an axis-based res_pair keyer_part with a cell_width that isn't a sensible, strictly-positive value"));
					}
				}

				/// \brief Get a short name that describes this key part
				static std::string get_name() {
					return "view_x";
				}

				/// \brief Extract the relevant value from the specified res_pair
				static view_base_type get_value(const multi_struc_res_rep_pair &arg_res_pair ///< The res_pair from which the relevant value should be extracted
				                                ) {
					return get_view_x( arg_res_pair.get_res_pair_core() );
				}

				/// \brief Extract the search radius from the specified quad_criteria
				static view_base_type get_search_radius(const quad_criteria &arg_criteria   ///< The criteria defining what is considered a match
				                                        ) {
					return std::sqrt( arg_criteria.get_maximum_squared_distance() );
				}
			};

		}

		/// \brief Type alias for keyer_part for x-axis of view
		using res_pair_view_x_keyer_part = detail::res_pair_view_axis_keyer_part<detail::res_pair_view_x_keyer_part_spec>;
	}
}

#endif
