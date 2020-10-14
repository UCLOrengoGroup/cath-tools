/// \file
/// \brief The get_sorting_scores class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_GET_SORTING_SCORES_HPP
#define _CATH_TOOLS_SOURCE_SRC_CLUSTAGGLOM_CLUSTAGGLOM_GET_SORTING_SCORES_HPP

#include "common/type_aliases.hpp"

namespace cath { namespace common { class id_of_str_bidirnl; } }

namespace cath {
	namespace clust {

		size_vec get_sorting_scores(const common::id_of_str_bidirnl &,
		                            const doub_vec &);

		size_vec get_sorting_scores(const common::id_of_str_bidirnl &);

	} // namespace clust
} // namespace cath

#endif
