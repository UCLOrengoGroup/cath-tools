/// \file
/// \brief The mapping_job class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_DETAIL_MAPPING_JOB_HPP
#define _CATH_TOOLS_SOURCE_CLUSTER_DETAIL_MAPPING_JOB_HPP

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "cluster/cluster_type_aliases.hpp"
#include "common/path_type_aliases.hpp"
#include "common/type_aliases.hpp"

namespace cath {
	namespace clust {
		namespace detail {

			/// \brief Represent a single mapping job to be performed
			///
			/// This is useful for doing batches of mapping (and in a common way)
			/// as for a single job
			class mapping_job final {
			private:
				/// \brief An optional identifier of the batch
				str_opt  batch_id;

				/// \brief The file describing the cluster membership of the clusters to be mapped/renumbered
				boost::filesystem::path new_cluster_membership_file;

				/// \brief An optional file describing map-from cluster membership
				path_opt old_cluster_membership_file;

			public:
				explicit mapping_job(const str_opt &,
				                     const boost::filesystem::path &,
				                     const path_opt & = boost::none);

				const str_opt & get_batch_id() const;
				const boost::filesystem::path & get_new_cluster_membership_file() const;
				const path_opt & get_old_cluster_membership_file() const;
			};

			mapping_job_vec read_batch_mapping_file(std::istream &);
			mapping_job_vec read_batch_mapping_file(const boost::filesystem::path &);

		} // namespace detail
	} // namespace clust
} // namespace cath

#endif
