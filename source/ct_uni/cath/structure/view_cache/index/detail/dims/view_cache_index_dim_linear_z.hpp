/// \file
/// \brief The view_cache_index_dim_linear class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_VIEW_CACHE_INDEX_DIM_LINEAR_Z_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_VIEW_CACHE_INDEX_DIM_LINEAR_Z_HPP

#include "cath/structure/view_cache/index/detail/dims/detail/vci_linear_dim_spec_view_axis.hpp"
#include "cath/structure/view_cache/index/detail/dims/detail/view_cache_index_dim_linear.hpp"

namespace cath::index::detail {

	namespace detail {

		/// \brief TODOCUMENT
		struct vci_linear_view_z_getter final {

			/// \brief TODOCUMENT
			[[nodiscard]] std::string get_name() const {
				return "z-axis";
			}

			/// \brief TODOCUMENT
			auto operator()(const view_cache_index_entry &prm_entry ///< TODOCUMENT
			                )->decltype( get_view_z( prm_entry ) ) {
				return get_view_z( prm_entry );
			}
		};

		using vci_linear_dim_spec_view_z = vci_linear_dim_spec_view_axis<vci_linear_view_z_getter>;

	} // namespace detail

	/// \brief TODOCUMENT
	using view_cache_index_dim_linear_z = detail::view_cache_index_dim_linear< detail::vci_linear_dim_spec_view_z >;

} // namespace cath::index::detail

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_VIEW_CACHE_INDEX_DIM_LINEAR_Z_HPP

