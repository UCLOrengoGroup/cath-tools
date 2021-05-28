/// \file
/// \brief The vci_linear_dim_spec_view_angle` class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL_VCI_LINEAR_DIM_SPEC_VIEW_ANGLE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL_VCI_LINEAR_DIM_SPEC_VIEW_ANGLE_HPP

#include "cath/common/config.hpp"
#include "cath/structure/view_cache/index/detail/dims/detail/vci_linear_dim_spec_view_angle.hpp"
#include "cath/structure/view_cache/index/detail/vcie_match_criteria.hpp"
#include "cath/structure/view_cache/index/detail/view_cache_index_type_aliases.hpp"

namespace cath::index::detail::detail {

	/// \brief TODOCUMENT
	template <typename F>
	struct vci_linear_dim_spec_view_angle final {
		/// \brief TODOCUMENT
		using value_type = angle_type;

		/// \brief TODOCUMENT
		using value_pair      = std::pair<value_type, value_type>;

		/// \brief TODOCUMENT
		using linear_dim_type = vci_linear_dim_spec_view_angle<vci_linear_dim_spec_view_angle>;

		/// \brief TODOCUMENT
		[[nodiscard]] std::string get_name() const {
			return F().get_name();
		}

		/// \brief TODOCUMENT
		void check_cell_width(const value_type &prm_cell_width ///< TODOCUMENT
		                      ) {
			if ( prm_cell_width <= geom::ZERO_ANGLE<angle_base_type> || prm_cell_width > geom::ONE_REVOLUTION<angle_base_type> ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Cannot create a vci_linear_dim_spec_view_angle<view_from_phi> with a cell_width that isn't in (0, 2pi ] "));
			}
		}

		/// \brief TODOCUMENT
		inline value_pair prepare_search_begin_and_end(value_type prm_search_min, ///< TODOCUMENT
		                                               value_type prm_search_max  ///< TODOCUMENT
		                                               ) {
			prm_search_min.quick_shift();
			prm_search_max.quick_shift();
			return std::make_pair( prm_search_min, prm_search_max );
		}

		/// \brief TODOCUMENT
		value_type get_search_radius(const vcie_match_criteria &prm_criteria ///< TODOCUMENT
		                             ) {
			if constexpr ( common::IS_IN_DEBUG_MODE ) {
				if ( F()( prm_criteria ) >= geom::HALF_REVOLUTION<angle_base_type> ) {
					BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Unable to search with radius >= pi"
						" (because then a wrapped end can't be reliably detected by checking whether it ends up less than the start)"));
				}
			}
			return F()( prm_criteria );
		}

		/// \brief TODOCUMENT
		const value_type & get_index_value(const view_cache_index_entry &prm_entry ///< TODOCUMENT
		                                   ) {
			const value_type &result = F()( prm_entry );
			if constexpr ( common::IS_IN_DEBUG_MODE ) {
				if ( result < geom::ZERO_ANGLE<angle_base_type> || result >= geom::ONE_REVOLUTION<angle_base_type> ) {
					BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Unable to index entry with angle out of range [0, 2pi)"));
				}
			}
			return result;
		}
	};

} // namespace cath::index::detail::detail

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL_VCI_LINEAR_DIM_SPEC_VIEW_ANGLE_HPP
