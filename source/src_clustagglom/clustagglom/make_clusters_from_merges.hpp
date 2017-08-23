/// \file
/// \brief The make_clusters_from_merges class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_AGGLOM_MAKE_CLUSTERS_FROM_MERGES_H
#define _CATH_TOOLS_SOURCE_CLUSTER_AGGLOM_MAKE_CLUSTERS_FROM_MERGES_H

#include <boost/filesystem/path.hpp>

#include "clustagglom/clustagglom_type_aliases.hpp"
#include "common/type_aliases.hpp"

namespace cath { namespace clust { struct hierarchy; } }
namespace cath { namespace common { class id_of_str_bidirnl; } }

namespace cath {
	namespace clust {
		namespace detail {

			merge_vec_citr_vec calc_merge_cutoff_boundaries(const merge_vec &,
			                                                const strength_vec &);

		} // namespace detail

		hierarchy make_clusters_from_merges(const merge_vec &,
		                                    const item_idx &,
		                                    const strength_vec &);

		hierarchy make_clusters_from_merges_and_sort(const merge_vec &,
		                                             const size_vec &,
		                                             const strength_vec &);


		void write_cluster(const boost::filesystem::path &,
		                   const hierarchy &,
		                   const common::id_of_str_bidirnl &);


	} // namespace clust
} // namespace cath

#endif
