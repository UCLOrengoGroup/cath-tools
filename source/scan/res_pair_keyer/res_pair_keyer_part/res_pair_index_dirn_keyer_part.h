/// \file
/// \brief The res_pair_index_dirn_keyer_part class header

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

#ifndef RES_PAIR_INDEX_DIRN_KEYER_PART_H_INCLUDED
#define RES_PAIR_INDEX_DIRN_KEYER_PART_H_INCLUDED

#include "boost_1_56_0/core/ignore_unused.hpp"

#include "exception/not_implemented_exception.h"

namespace cath {
	namespace scan {

		/// \brief Key part generator for multi_struc_res_rep_pair that extracts the index direction
		class res_pair_index_dirn_keyer_part final {
		public:
			/// \brief TODOCUMENT
			using value_t           = detail::res_pair_dirn;

			/// \brief TODOCUMENT
			using cell_index_t      = detail::res_pair_dirn;

			/// \brief TODOCUMENT
			using cell_index_list_t = std::array<detail::res_pair_dirn, 1>;

			/// \brief TODOCUMENT
			using search_radius_t   = res_pair_index_dirn_criterion;

			/// \brief Get a short name that describes this key part
			std::string get_name() const {
				return "index_dirn";
			}

			/// \brief Extract the relevant value from the specified res_pair
			value_t get_value(const detail::multi_struc_res_rep_pair &arg_res_pair ///< The res_pair from which the relevant value should be extracted
			                  ) const {
				return direction( arg_res_pair );
			}

			/// \brief Extract the search radius from the specified quad_criteria
			search_radius_t get_search_radius(const quad_criteria &arg_criteria ///< The criteria defining what is considered a match
			                                  ) const {
				return arg_criteria.get_index_direction_criterion();
			}

			/// \brief Generate the key part for the specified value
			cell_index_t key_part(const value_t &arg_value ///< The value for which the key_part should be extracted
			                      ) const {
				return arg_value;
			}

			/// \brief Generate a list of all key parts for all conceivable res_pairs that would match the specified value
			///        within the specified search radius
			cell_index_list_t close_key_parts(const value_t         &arg_value,        ///< The value for which the key_part should be extracted
			                                  const search_radius_t &arg_search_radius ///< The search radius defining what is considered a match
			                                  ) const {
#ifndef NDEBUG
				if ( arg_search_radius != res_pair_index_dirn_criterion::MUST_MATCH ) {
					BOOST_THROW_EXCEPTION(common::not_implemented_exception("res_pair_index_dirn_keyer_part currently requires that the search radius is false (ie require_matching_directions)"));
				}
#else
				boost::ignore_unused( arg_search_radius );
#endif
				return { { arg_value } };
			}
		};

	}
}

#endif
