/// \file
/// \brief The res_pair_from_to_index_keyer_part class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_FROM_TO_INDEX_KEYER_PART_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_FROM_TO_INDEX_KEYER_PART_HPP

#include "cath/common/config.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/common/tuple/make_tuple_with_skips.hpp"
#include "cath/scan/spatial_index/simple_locn_index.hpp"

namespace cath::scan {

	/// \brief Key part generator for multi_struc_res_rep_pair that extracts the index direction
	class res_pair_from_to_index_keyer_part final {
	public:
		/// \brief TODOCUMENT
		using value_t           = unsigned int;

		/// \brief TODOCUMENT
		using cell_index_t      = common::tpl_elmnt_skip_t;

		/// \brief TODOCUMENT
		using cell_index_list_t = common::tpl_elmnt_skip_t;

		/// \brief TODOCUMENT
		using search_radius_t   = unsigned int;

		/// \brief Get a short name that describes this key part
		[[nodiscard]] std::string get_name() const {
			return "from_to_index";
		}

		/// \brief Extract the relevant value from the specified res_pair
		[[nodiscard]] constexpr value_t get_value(const simple_locn_index &prm_res_pair ///< The res_pair from which the relevant value should be extracted
		                                          ) const {
			return prm_res_pair.index;
		}

		/// \brief Extract the search radius from the specified simple_locn_crit
		[[nodiscard]] constexpr search_radius_t get_search_radius(const simple_locn_crit &/*prm_criteria*/ ///< The criteria defining what is considered a match
		                                                          ) const {
			return 0;
		}

		/// \brief Generate the key part for the specified value
		[[nodiscard]] constexpr cell_index_t key_part(const value_t &/*prm_value*/ ///< The value for which the key_part should be extracted
		                                              ) const {
			return {};
		}

		/// \brief Generate a list of all key parts for all conceivable res_pairs that would match the specified value
		///        within the specified search radius
		[[nodiscard]] cell_index_list_t close_key_parts(const value_t         &/*prm_value*/,        ///< The value for which the key_part should be extracted
		                                                const search_radius_t &prm_search_radius ///< The search radius defining what is considered a match
		                                                ) const {
			if constexpr ( common::IS_IN_DEBUG_MODE ) {
				if ( prm_search_radius != 0 ) {
					BOOST_THROW_EXCEPTION( common::not_implemented_exception(
					  "res_pair_from_to_index_keyer_part currently requires that the search radius is 0" ) );
				}
			}
			return {};
		}

		/// \brief Generate the minimum key part within the specified search radius for the specified value
		///
		/// This is an extra that makes this usable for dense storing in a lattice
		[[nodiscard]] constexpr cell_index_t min_close_key_part(const value_t         &/*prm_value*/,    ///< The value for which the key_part should be extracted
		                                                        const search_radius_t &prm_search_radius ///< The search radius defining what is considered a match
		                                                        ) const {
			if constexpr ( common::IS_IN_DEBUG_MODE ) {
				if ( prm_search_radius != 0 ) {
					BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
					  "res_pair_from_to_index_keyer_part currently requires that the search radius is 0" ) );
				}
			}
			return cell_index_t{};
		}

		/// \brief Generate the maximum key part within the specified search radius for the specified value
		///
		/// This is an extra that makes this usable for dense storing in a lattice
		[[nodiscard]] constexpr cell_index_t max_close_key_part(const value_t         &/*prm_value*/,    ///< The value for which the key_part should be extracted
		                                                        const search_radius_t &prm_search_radius ///< The search radius defining what is considered a match
		                                                        ) const {
			if constexpr ( common::IS_IN_DEBUG_MODE ) {
				if ( prm_search_radius != 0 ) {
					BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
					  "res_pair_from_to_index_keyer_part currently requires that the search radius is 0" ) );
				}
			}
			return cell_index_t{};
		}
	};

} // namespace cath::scan

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_RES_PAIR_KEYER_RES_PAIR_KEYER_PART_RES_PAIR_FROM_TO_INDEX_KEYER_PART_HPP
