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

#ifndef _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_CALC_COMPLETE_LINKAGE_MERGE_LIST_H
#define _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_CALC_COMPLETE_LINKAGE_MERGE_LIST_H

#include "clustagglom/clustagglom_type_aliases.hpp"
#include "common/type_aliases.hpp"

namespace cath { namespace clust { class links; } }

namespace cath {
	namespace clust {

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

	} // namespace clust
} // namespace cath

#endif
