/// \file
/// \brief The cath_cluster_clustering_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC_CATH_CLUSTER_CLUSTERING_SPEC_HPP
#define _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC_CATH_CLUSTER_CLUSTERING_SPEC_HPP

#include "cath/clustagglom/clustagglom_type_aliases.hpp"
#include "cath/clustagglom/link_dirn.hpp"
#include "cath/common/path_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath::clust {

	/// \brief Specify the output for cath-cluster
	class cath_cluster_clustering_spec final {
	private:
		/// \brief The levels at which the clustering should be performed
		strength_vec levels;

	  public:
		[[nodiscard]] const strength_vec &get_levels() const;

		cath_cluster_clustering_spec &set_levels( const strength_vec & );
	};

	str_opt get_invalid_description(const cath_cluster_clustering_spec &);



	str_opt get_dissim_sort_warning(const strength_vec &,
	                                const link_dirn &);

	void make_dissim_and_sort(strength_vec &,
	                          const link_dirn &);

	strength_vec make_dissim_and_sort_copy(strength_vec,
	                                       const link_dirn &);


	strength get_max_dissim(const strength_vec &,
	                        const link_dirn &);

	strength_vec get_sorted_dissims(const cath_cluster_clustering_spec &,
	                                const link_dirn &);

	str_opt get_dissim_sort_warning(const cath_cluster_clustering_spec &,
	                                const link_dirn &);


	strength get_max_dissim(const cath_cluster_clustering_spec &,
	                        const link_dirn &);

} // namespace cath::clust

#endif // _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_SPEC_CATH_CLUSTER_CLUSTERING_SPEC_HPP
