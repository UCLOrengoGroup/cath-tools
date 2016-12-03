/// \file
/// \brief The axis_keyer_part class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_DETAIL_AXIS_KEYER_PART_H
#define _CATH_TOOLS_SOURCE_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_DETAIL_AXIS_KEYER_PART_H

#include <boost/range/irange.hpp>

#include "common/algorithm/constexpr_floor.hpp"
#include "common/debug_numeric_cast.hpp"
#include "scan/detail/scan_type_aliases.hpp"

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief Template for keyer_part that indexes the view's axes
			template <typename Spec>
			class axis_keyer_part final {
			public:
				/// \brief TODOCUMENT
				using value_t           = typename Spec::value_t;

				/// \brief TODOCUMENT
				using stored_t          = typename Spec::stored_t;

				/// \brief TODOCUMENT
				using criteria_t        = typename Spec::criteria_t;

				/// \brief TODOCUMENT
				using cell_index_t      = key_view_index_type;

				/// \brief TODOCUMENT
				using cell_index_list_t = boost::integer_range<cell_index_t>;

				/// \brief TODOCUMENT
				using search_radius_t   = value_t;

			private:
				/// \brief The cell width to use for this dimension
				value_t cell_width;

			public:
				/// \brief Ctor from cell_width
				explicit constexpr axis_keyer_part(const value_t &arg_cell_width ///< The cell width to use in keying this part
				                                   ) : cell_width( Spec::sanity_check_cell_width( arg_cell_width ) ) {
				}

				/// \brief Get a short name that describes this key part
				std::string get_name() const {
					return Spec::get_name() + "(" + std::to_string( cell_width ) + "]";
				}

				/// \brief Extract the relevant value from the specified res_pair
				constexpr value_t get_value(const stored_t &arg_res_pair ///< The res_pair from which the relevant value should be extracted
				                            ) const {
					return Spec::get_value( arg_res_pair );
				}

				/// \brief Extract the search radius from the specified criteria
				constexpr search_radius_t get_search_radius(const criteria_t &arg_criteria   ///< The criteria defining what is considered a match
				                                            ) const {
					return Spec::get_search_radius( arg_criteria  );
				}

				/// \brief Generate the key part for the specified value
				constexpr cell_index_t key_part(const value_t &arg_value ///< The value for which the key_part should be extracted
				                                ) const {
					return static_cast<cell_index_t>( common::constexpr_floor( arg_value / cell_width ) );
				}

				/// \brief Generate a list of all key parts for all conceivable res_pairs that would match the specified value
				///        within the specified search radius
				cell_index_list_t close_key_parts(const value_t         &arg_value,        ///< The res_pair whose matches' key parts should be generated
				                                  const search_radius_t &arg_search_radius ///< The search radius defining what is considered a match
				                                  ) const {
					return boost::irange(
						key_part( arg_value - arg_search_radius ),
						debug_numeric_cast<cell_index_t>( key_part( arg_value + arg_search_radius ) + 1 )
					);
				}
			};


		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
