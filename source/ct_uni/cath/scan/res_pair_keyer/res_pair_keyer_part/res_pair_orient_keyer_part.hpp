/// \file
/// \brief The res_pair_orient_keyer_part class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_ORIENT_KEYER_PART_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_ORIENT_KEYER_PART_HPP

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

#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/scan/detail/res_pair/multi_struc_res_rep_pair.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"
#include "cath/structure/geometry/orientation_covering.hpp"

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Key part generator for multi_struc_res_rep_pair that indexes the view's x-dimension
			///
			/// !!!!!!!!! THIS CLASS IS ALMOST CERTAINLY HALF-WRITTEN !!!!!!!!!
			/// !!!!!!!!! THIS CLASS IS ALMOST CERTAINLY HALF-WRITTEN !!!!!!!!!
			/// !!!!!!!!! THIS CLASS IS ALMOST CERTAINLY HALF-WRITTEN !!!!!!!!!
			class res_pair_orient_keyer_part final {
			private:
				// /// \brief TODOCUMENT
				// geom::orientation_covering_impl<quat_rot_type> the_orientation_covering;

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
				[[nodiscard]] key_view_index_type get_key_part_of_value(const view_base_type &prm_value ///< The value for which the key part should be calculated
				                                          ) const {
					BOOST_THROW_EXCEPTION(common::not_implemented_exception("!!!!!!!!! THIS CLASS IS ALMOST CERTAINLY HALF-WRITTEN !!!!!!!!!"));
					return static_cast<key_view_index_type>( floor( prm_value / cell_width ) );
				}

				/// \brief Extract the relevant value from the specified res_pair
				[[nodiscard]] view_base_type get_value(const multi_struc_res_rep_pair &prm_res_pair ///< The res_pair from which the relevant value should be extracted
				                         ) const {
					BOOST_THROW_EXCEPTION(common::not_implemented_exception("!!!!!!!!! THIS CLASS IS ALMOST CERTAINLY HALF-WRITTEN !!!!!!!!!"));
					return get_view_x( prm_res_pair.get_res_pair_core() );
				}

			public:
				/// \brief Generate the key part for the specified res_pair
				[[nodiscard]] key_view_index_type key_part(const view_base_type &prm_value ///< The value for which the key_part should be extracted
				                             // const search_radius_t &prm_search_radius ///< The search radius defining what is considered a match
				                             ) const {
					BOOST_THROW_EXCEPTION(common::not_implemented_exception("!!!!!!!!! THIS CLASS IS ALMOST CERTAINLY HALF-WRITTEN !!!!!!!!!"));
					return get_key_part_of_value( prm_value );
				}

				/// \brief Generate a list of all key parts for all conceivable res_pairs that would match the specified res_pair
				[[nodiscard]] boost::integer_range<key_view_index_type> close_key_parts(const multi_struc_res_rep_pair &prm_res_pair,     ///< The res_pair whose matches' key parts should be generated
				                                                          const search_radius_t          &prm_search_radius ///< The criteria defining what is considered a match
				                                                          ) const {
					BOOST_THROW_EXCEPTION(common::not_implemented_exception("!!!!!!!!! THIS CLASS IS ALMOST CERTAINLY HALF-WRITTEN !!!!!!!!!"));
					return boost::irange<key_view_index_type>(
						                                         get_key_part_of_value( get_value( prm_res_pair ) - prm_search_radius ),
						debug_numeric_cast<key_view_index_type>( get_key_part_of_value( get_value( prm_res_pair ) + prm_search_radius ) + 1 )
					);
				}
			};

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_ORIENT_KEYER_PART_HPP
