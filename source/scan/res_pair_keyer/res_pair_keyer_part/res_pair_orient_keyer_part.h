/// \file
/// \brief The res_pair_orient_keyer_part class header

#ifndef RES_PAIR_ORIENT_KEYER_PART_H_INCLUDED
#define RES_PAIR_ORIENT_KEYER_PART_H_INCLUDED

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "structure/geometry/orientation_covering.h"

namespace cath {
	namespace scan {

		/// \brief Key part generator for multi_struc_res_rep_pair that indexes the view's x-dimension
		class res_pair_orient_keyer_part final {
		private:
			/// \brief TODOCUMENT
			geom::orientation_covering the_orientation_covering;

			/// \brief TODOCUMENT
			using value_t           = frame_quat_rot;

			/// \brief TODOCUMENT
			using cell_index_t      = key_orient_index_type;

			/// \brief TODOCUMENT
			using cell_index_list_t = std::vector<key_orient_index_type>;

			/// \brief TODOCUMENT
			using search_radius_t   = frame_quat_rot_type;

			/// \brief The cell width to use for the view's x-dimension
			view_base_type cell_width;

			/// \brief Calculate the key part (ie cell index) for the specified value
			key_view_index_type get_key_part_of_value(const view_base_type &arg_value ///< The value for which the key part should be calculated
			                                          ) const {
				return static_cast<key_view_index_type>( floor( value / cell_width ) );
			}

			/// \brief Extract the relevant value from the specified res_pair
			view_base_type get_value(const multi_struc_res_rep_pair &arg_res_pair ///< The res_pair from which the relevant value should be extracted
			                         ) const {
				return get_view_x( arg_res_pair.get_res_pair_core() );
			}

		public:
			/// \brief Generate the key part for the specified res_pair
			key_view_index_type key_part(const value_t         &arg_value,        ///< The value for which the key_part should be extracted
			                             const search_radius_t &arg_search_radius ///< The search radius defining what is considered a match
			                             ) const {
				return get_key_part_of_value( get_value( arg_res_pair ) );
			}

			/// \brief Generate a list of all key parts for all conceivable res_pairs that would match the specified res_pair
			boost::integer_range<key_view_index_type> close_key_parts(const multi_struc_res_rep_pair &arg_res_pair, ///< The res_pair whose matches' key parts should be generated
			                                                          const search_radius_t          &arg_search_radius ///< The criteria defining what is considered a match
			                                                          ) const {
				const auto search_radius = std::sqrt( arg_criteria.get_maximum_squared_distance() );
				return irange(
					get_key_part_of_value( get_value( arg_res_pair ) - search_radius ),
					get_key_part_of_value( get_value( arg_res_pair ) + search_radius ) + 1
				);
			}
		};


	}
}

#endif
