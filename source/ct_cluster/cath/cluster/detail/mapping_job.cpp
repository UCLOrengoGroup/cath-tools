/// \file
/// \brief The mapping_job class definitions

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

#include "mapping_job.hpp"

#include <boost/algorithm/string/classification.hpp>

#include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/optional/make_optional_if.hpp"

#include <fstream>

using namespace ::cath;
using namespace ::cath::clust::detail;
using namespace ::cath::common;

using ::boost::algorithm::is_any_of;
using ::boost::filesystem::path;
using ::boost::token_compress_on;
using ::std::ifstream;
using ::std::istream;
using ::std::string;

/// \brief Ctor
mapping_job::mapping_job(const str_opt  &prm_batch_id,                    ///< An optional identifier of the batch
                         const path     &prm_new_cluster_membership_file, ///< The file describing the cluster membership of the clusters to be mapped/renumbered
                         const path_opt &prm_old_cluster_membership_file  ///< An optional file describing map-from cluster membership
                         ) : batch_id                   { prm_batch_id                    },
                             new_cluster_membership_file{ prm_new_cluster_membership_file },
                             old_cluster_membership_file{ prm_old_cluster_membership_file } {
}

/// \brief Getter for an optional identifier of the batch
const str_opt & mapping_job::get_batch_id() const {
	return batch_id;
}

/// \brief Getter for the file describing the cluster membership of the clusters to be mapped/renumbered
const path & mapping_job::get_new_cluster_membership_file() const {
	return new_cluster_membership_file;
}

/// \brief Getter for an optional file describing map-from cluster membership
const path_opt & mapping_job::get_old_cluster_membership_file() const {
	return old_cluster_membership_file;
}

/// \brief Read a batch of mapping jobs from the specified istream
///
/// \brief relates mapping_job
mapping_job_vec cath::clust::detail::read_batch_mapping_file(istream &prm_istream ///< The istream from which to read the mapping jobs
                                                             ) {
	mapping_job_vec results;
	string line;
	while ( getline( prm_istream, line ) ) {
		const auto parts = split_build<str_vec>( line, is_any_of( " \t" ), token_compress_on );
		if ( parts.size() < 2 ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("There should be at least two columns (batch_id working_clust_memb_file) in batch specification input lines"));
		}
		if ( parts.size() > 3 ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("There cannot be more than three columns (batch_id working_clust_memb_file prev_clust_memb_file) in batch specification input lines"));
		}
		string batch_id_str( parts.front() );
		if ( batch_id_str == "-" ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Cannot use files with name \"-\" in batch specification input lines"));
		}

		results.emplace_back(
			move( batch_id_str ),
			path{ parts[ 1 ] },
			make_optional_if_fn(
				parts.size() >= 3,
				[&] { return path{ parts[ 2 ] }; }
			)
		);
	}
	return results;
}

/// \brief Read a batch of mapping jobs from the specified file
///
/// \brief relates mapping_job
mapping_job_vec cath::clust::detail::read_batch_mapping_file(const path &prm_file ///< The file from which to read the mapping jobs
                                                             ) {
	ifstream mapping_istream;
	open_ifstream( mapping_istream, prm_file );
	mapping_job_vec results = read_batch_mapping_file( mapping_istream );
	mapping_istream.close();
	return results;
}
