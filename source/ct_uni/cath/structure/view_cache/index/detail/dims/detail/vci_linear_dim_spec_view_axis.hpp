/// \file
/// \brief The vci_linear_dim_spec_view_axis class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL_VCI_LINEAR_DIM_SPEC_VIEW_AXIS_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL_VCI_LINEAR_DIM_SPEC_VIEW_AXIS_HPP

#include "cath/structure/view_cache/index/detail/dims/detail/view_cache_index_dim_linear.hpp"
#include "cath/structure/view_cache/index/detail/view_cache_index_type_aliases.hpp"

namespace cath::index::detail::detail {

	/// \brief TODOCUMENT
	template <typename F>
	struct vci_linear_dim_spec_view_axis final {
		/// \brief TODOCUMENT
		using value_type      = view_base_type;

		/// \brief TODOCUMENT
		using value_pair      = std::pair<value_type, value_type>;

		/// \brief TODOCUMENT
		using linear_dim_type = view_cache_index_dim_linear<vci_linear_dim_spec_view_axis>;

		/// \brief TODOCUMENT
		[[nodiscard]] std::string get_name() const {
			return F().get_name();
		}

		/// \brief TODOCUMENT
		void check_cell_width(const value_type &prm_cell_width ///< TODOCUMENT
		                      ) {
			if ( ! boost::math::isnormal( prm_cell_width ) || prm_cell_width <= 0.0 ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Cannot create a view_cache_index_dim_linear<vci_linear_dim_spec_view_axis<>> with a cell_width that isn't a sensible, strictly-positive value"));
			}
		}

		/// \brief TODOCUMENT
		inline value_pair prepare_search_begin_and_end(const value_type &prm_search_min, ///< TODOCUMENT
		                                               const value_type &prm_search_max  ///< TODOCUMENT
		                                               ) {
			// std::cerr << "In vci_linear_dim_spec_view_axis::prepare_search_begin_and_end(), prm_search_min is " << prm_search_min << std::endl;
			// std::cerr << "In vci_linear_dim_spec_view_axis::prepare_search_begin_and_end(), prm_search_max is " << prm_search_max << std::endl;
			// return std::make_tuple(
			// 	true,
			// 	prm_linear_dim.cell_index_of_value_in_current( prm_search_min ),
			// 	prm_linear_dim.cell_index_of_value_in_current( prm_search_max ) + 1
			// );
			assert( prm_search_min < prm_search_max );
			return std::make_pair( prm_search_min, prm_search_max );
		}

		/// \brief TODOCUMENT
		value_type get_search_radius(const vcie_match_criteria &prm_criteria ///< TODOCUMENT
		                             ) {
			return std::sqrt( prm_criteria.get_maximum_squared_distance() );
		}

		/// \brief TODOCUMENT
		const value_type & get_index_value(const view_cache_index_entry &prm_entry ///< TODOCUMENT
		                                   ) {
			return F()( prm_entry );
		}
	};

} // namespace cath::index::detail::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL_VCI_LINEAR_DIM_SPEC_VIEW_AXIS_HPP

