/// \file
/// \brief The res_pair_index_dirn_keyer_part class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_INDEX_DIRN_KEYER_PART_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_INDEX_DIRN_KEYER_PART_HPP

#include "cath/common/config.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/scan/detail/res_pair/multi_struc_res_rep_pair.hpp"
#include "cath/scan/detail/res_pair_dirn/res_pair_dirn.hpp"
#include "cath/scan/quad_criteria.hpp"
#include "cath/scan/res_pair_index_dirn_criterion.hpp"

namespace cath::scan {

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
		[[nodiscard]] std::string get_name() const {
			return "index_dirn";
		}

		/// \brief Extract the relevant value from the specified res_pair
		[[nodiscard]] value_t get_value(const detail::multi_struc_res_rep_pair &prm_res_pair ///< The res_pair from which the relevant value should be extracted
		                  ) const {
			return direction( prm_res_pair );
		}

		/// \brief Extract the search radius from the specified quad_criteria
		[[nodiscard]] search_radius_t get_search_radius(const quad_criteria &prm_criteria ///< The criteria defining what is considered a match
		                                  ) const {
			return prm_criteria.get_index_direction_criterion();
		}

		/// \brief Generate the key part for the specified value
		[[nodiscard]] cell_index_t key_part(const value_t &prm_value ///< The value for which the key_part should be extracted
		                      ) const {
			return prm_value;
		}

		/// \brief Generate a list of all key parts for all conceivable res_pairs that would match the specified value
		///        within the specified search radius
		[[nodiscard]] cell_index_list_t close_key_parts(const value_t         &prm_value,        ///< The value for which the key_part should be extracted
		                                  const search_radius_t &prm_search_radius ///< The search radius defining what is considered a match
		                                  ) const {
			if constexpr ( common::IS_IN_DEBUG_MODE ) {
				if ( prm_search_radius != res_pair_index_dirn_criterion::MUST_MATCH ) {
					BOOST_THROW_EXCEPTION(
					  common::not_implemented_exception( "res_pair_index_dirn_keyer_part currently requires that the "
					                                     "search radius is false (ie require_matching_directions)" ) );
				}
			}
			return { { prm_value } };
		}
	};

} // namespace cath::scan

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_INDEX_DIRN_KEYER_PART_HPP
