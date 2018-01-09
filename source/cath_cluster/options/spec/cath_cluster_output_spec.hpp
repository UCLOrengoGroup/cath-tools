/// \file
/// \brief The cath_cluster_output_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_CATH_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CATH_CLUSTER_OUTPUT_SPEC_HPP
#define _CATH_TOOLS_SOURCE_SRC_CATH_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CATH_CLUSTER_OUTPUT_SPEC_HPP

#include <boost/optional.hpp>

#include "common/path_type_aliases.hpp"
#include "common/type_aliases.hpp"

namespace cath {
	namespace clust {

		/// \brief Specify the output for cath-cluster
		class cath_cluster_output_spec final {
		private:
			/// \brief An optional file to which clusters should be written
			path_opt clusters_to_file;

			/// \brief An optional file to which merges should be written
			path_opt merges_to_file;

			/// \brief An optional file to which clust_spans should be written
			path_opt clust_spans_to_file;

			/// \brief An optional file to which reps should be written
			path_opt reps_to_file;

			/// \brief An optional file to which sorted_links should be written
			path_opt sorted_links_to_file;

		public:
			const path_opt & get_clusters_to_file() const;
			const path_opt & get_merges_to_file() const;
			const path_opt & get_clust_spans_to_file() const;
			const path_opt & get_reps_to_file() const;
			const path_opt & get_sorted_links_to_file() const;

			cath_cluster_output_spec & set_clusters_to_file(const path_opt &);
			cath_cluster_output_spec & set_merges_to_file(const path_opt &);
			cath_cluster_output_spec & set_clust_spans_to_file(const path_opt &);
			cath_cluster_output_spec & set_reps_to_file(const path_opt &);
			cath_cluster_output_spec & set_sorted_links_to_file(const path_opt &);
		};

		path_vec get_all_output_paths(const cath_cluster_output_spec &);

		str_opt get_invalid_description(const cath_cluster_output_spec &);

	} // namespace clust
} // namespace cath

#endif
