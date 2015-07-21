/// \file
/// \brief The view_cache_index_dim_linear_from_phi class header

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

#ifndef VIEW_CACHE_INDEX_DIM_LINEAR_FROM_PHI_H_INCLUDED
#define VIEW_CACHE_INDEX_DIM_LINEAR_FROM_PHI_H_INCLUDED

#include "structure/view_cache/index/detail/dims/detail/vci_linear_dim_spec_view_angle.h"
#include "structure/view_cache/index/detail/dims/detail/view_cache_index_dim_linear.h"

namespace cath {
	namespace index {
		namespace detail {
			namespace detail {
				
				/// \brief TODOCUMENT
				struct vci_linear_view_from_phi_getter final {

					/// \brief TODOCUMENT
					std::string get_name() const {
						return "from-phi";
					}

					/// \brief TODOCUMENT
					auto operator()(const view_cache_index_entry &arg_entry ///< TODOCUMENT
					                )->decltype( arg_entry.get_from_phi_angle() ) {
						return arg_entry.get_from_phi_angle();
					}

					/// \brief TODOCUMENT
					auto operator()(const vcie_match_criteria &arg_criteria ///< TODOCUMENT
					                )->decltype( arg_criteria.get_maximum_phi_angle_difference() ) {
						return arg_criteria.get_maximum_phi_angle_difference();
					}
				};

				using vci_linear_dim_spec_view_from_phi = vci_linear_dim_spec_view_angle<vci_linear_view_from_phi_getter>;
			}

			/// \brief TODOCUMENT
			using view_cache_index_dim_linear_from_phi = detail::view_cache_index_dim_linear<detail::vci_linear_dim_spec_view_from_phi>;
		}
	}
}

#endif

