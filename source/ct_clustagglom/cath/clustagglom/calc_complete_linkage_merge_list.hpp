/// \file
/// \brief The calc_complete_linkage_merge_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_CALC_COMPLETE_LINKAGE_MERGE_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_CALC_COMPLETE_LINKAGE_MERGE_LIST_HPP

#include "cath/clustagglom/clustagglom_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::clust { class links; }
// clang-format on

namespace cath::clust {

	merge_vec calc_complete_linkage_merge_list(links,
	                                           const size_vec &,
	                                           const strength & = std::numeric_limits<strength>::infinity() );

	merge_vec calc_complete_linkage_merge_list(links,
	                                           const size_t &,
	                                           const strength & = std::numeric_limits<strength>::infinity() );

	merge_vec calc_complete_linkage_merge_list(const item_item_strength_tpl_vec &,
	                                           const size_vec &,
	                                           const strength & = std::numeric_limits<strength>::infinity() );

	merge_vec calc_complete_linkage_merge_list(const item_item_strength_tpl_vec &,
	                                           const size_t &,
	                                           const strength & = std::numeric_limits<strength>::infinity() );

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_CALC_COMPLETE_LINKAGE_MERGE_LIST_HPP
