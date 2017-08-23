/// \map
/// \brief The map_results class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_MAP_MAP_RESULTS_H
#define _CATH_TOOLS_SOURCE_CLUSTER_MAP_MAP_RESULTS_H

#include <boost/filesystem/path.hpp>

#include "cluster/cluster_type_aliases.hpp"
#include "cluster/map/overlap_frac_distn.hpp"
#include "cluster/options/spec/clust_mapping_spec.hpp"
#include "common/type_aliases.hpp"

#include <iostream>
#include <string>

namespace cath { namespace clust { class clust_mapping_spec; } }
namespace cath { namespace clust { class new_cluster_data; } }

namespace cath {
	namespace clust {

		/// \brief Store the basic data on a potential map between to clusters
		struct potential_map final {
			/// \brief The index of the old, map-from cluster
			size_t old_cluster_idx;

			/// \brief The index of the new, map-to cluster
			size_t new_cluster_idx;

			/// \brief The number of entries that map between the two clusters
			size_t num_mapped;

			/// \brief Natural ctor to populate each element
			potential_map(const size_t &arg_old_cluster_idx, ///< The index of the old, map-from cluster
			              const size_t &arg_new_cluster_idx, ///< The index of the new, map-to cluster
			              const size_t &arg_num_mapped       ///< The number of entries that map between the two clusters
			              ) : old_cluster_idx { arg_old_cluster_idx },
			                  new_cluster_idx { arg_new_cluster_idx },
			                  num_mapped      { arg_num_mapped      } {
			}
		};

		namespace detail {
			std::string get_name_of_new_unmapped_cluster_of_index(const boost::optional<ptrdiff_t> &,
			                                                      const size_t &);

		}

		/// \brief Store the data arising from a mapping
		struct map_results final {
			/// \brief The maps that have finally been selected
			///
			/// These should mention each cluster at most once
			potential_map_vec chosen_maps;

			/// \brief The remaining maps that weren't chosen
			potential_map_vec other_maps;

			/// \brief The indices of the unmapped, new clusters
			///        sorted according to the cluster_info criteria
			size_vec unmapped_new_cluster_indices;

			/// \brief The total number of mapped entries for each of the new clusters, indexed by their clusters' indices
			size_vec num_mapped_by_new_cluster;

			/// \brief The total number of mapped entries for each of the old clusters, indexed by their clusters' indices
			size_vec num_mapped_by_old_cluster;

			/// \brief The number of entries which got a domain overlap of 0 *because there were no matching entries on the parent at all*
			size_t num_with_nothing_on_parent;

			/// \brief The highest overlap fraction (over largest) for each of the old domains
			overlap_frac_distn highest_old_dom_overlap_fractions;

			/// \brief The highest overlap fraction for each of the old clusters
			overlap_frac_distn highest_old_clust_overlap_fractions;

			/// \brief The spec under which the mapping was performed
			clust_mapping_spec the_spec;
		};

		size_t get_num_mapped_entries(const map_results &);

		std::string results_string(const old_cluster_data_opt &,
		                           const new_cluster_data &,
		                           const map_results &,
		                           const str_opt &,
		                           const bool & = true);

		std::string longer_results_string(const old_cluster_data_opt &,
		                                  const new_cluster_data &,
		                                  const map_results &,
		                                  const str_opt &);

		std::string markdown_summary_string(const old_cluster_data_opt &,
		                                    const new_cluster_data &,
		                                    const map_results &);

		void write_markdown_summary_string_to_file(const boost::filesystem::path &,
		                                           const old_cluster_data_opt &,
		                                           const new_cluster_data &,
		                                           const map_results &);

	} // namespace clust
} // namespace cath

#endif
