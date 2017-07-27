/// \map
/// \brief The aggregate_map_results class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_MAP_AGGREGATE_MAP_RESULTS_H
#define _CATH_TOOLS_SOURCE_CLUSTER_MAP_AGGREGATE_MAP_RESULTS_H

#include <boost/filesystem/path.hpp>

#include "cluster/map/overlap_frac_distn.hpp"
#include "cluster/options/spec/clust_mapping_spec.hpp"

namespace cath { namespace clust { class new_cluster_data; } }
namespace cath { namespace clust { class old_cluster_data; } }
namespace cath { namespace clust { struct map_results; } }

namespace cath {
	namespace clust {

		/// \brief Store the data summarising the results of one or more mappings
		///
		/// This requires that all mapping were performed with the same clust_mapping_spec
		class aggregate_map_results final {
		private:
			bool               added_to                 = false;

			/// \brief The total number of old clusters encountered in the mapping
			size_t             num_old_clusters         = 0;

			/// \brief The total number of new clusters encountered in the mapping
			size_t             num_new_clusters         = 0;

			/// \brief The total number of clusters (on either side) that got mapped
			///        (which is the same as the number of cluster mappings because each cluster can be mapped at most once)
			size_t             num_mapped_clusters      = 0;

			/// \brief The total number of old entries encountered in the mapping
			size_t             num_old_entries          = 0;

			/// \brief The total number of new entries encountered in the mapping
			size_t             num_new_entries          = 0;

			/// \brief The total number of entries on either side that got mapped
			///        (which is the same as the number of entry mappings because each entry can be mapped at most once)
			size_t             num_mapped_entries       = 0;

			/// \brief The highest overlap fraction (over largest) for each of the old domains
			overlap_frac_distn highest_old_dom_overlap_fractions;

			/// \brief The highest overlap fraction for each of the old clusters
			overlap_frac_distn highest_old_clust_overlap_fractions;

			/// \brief The specification that was used to perform the mappings
			clust_mapping_spec the_spec;

		public:
			aggregate_map_results() = default;
			explicit aggregate_map_results(const clust_mapping_spec &) noexcept;

			const bool & get_added_to() const;
			const size_t & get_num_old_clusters() const;
			const size_t & get_num_new_clusters() const;
			const size_t & get_num_mapped_clusters() const;
			const size_t & get_num_old_entries() const;
			const size_t & get_num_new_entries() const;
			const size_t & get_num_mapped_entries() const;
			const overlap_frac_distn & get_highest_old_dom_overlap_fractions() const;
			const overlap_frac_distn & get_highest_old_clust_overlap_fractions() const;
			const clust_mapping_spec & get_clust_mapping_spec() const;

			aggregate_map_results & add_map_results(const map_results &,
			                                        const old_cluster_data &,
			                                        const new_cluster_data &);
		};

		aggregate_map_results make_aggregate_map_results(const map_results &,
		                                                 const old_cluster_data &,
		                                                 const new_cluster_data &);

		std::string markdown_summary_string(const aggregate_map_results &);

		void write_markdown_summary_string_to_file(const boost::filesystem::path &,
		                                           const aggregate_map_results &);

	} // namespace clust
} // namespace cath

#endif
