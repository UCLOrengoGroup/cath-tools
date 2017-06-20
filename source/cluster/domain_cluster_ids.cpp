/// \file
/// \brief The domain_cluster_ids class definitions

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

#include "domain_cluster_ids.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <iostream>
#include <string>

using namespace cath::clust;

using boost::adaptors::transformed;
using boost::algorithm::join;
using std::string;

/// \brief Generate a string describing the specified domain_cluster_id
///
/// \relates domain_cluster_id
string cath::clust::to_string(const domain_cluster_id &arg_dom_clust_id,   ///< The domain_cluster_id to describe
                              const bool              &arg_include_cluster ///< Whether to include the cluster ID for each entry (default: true)
                              ) {
	using std::to_string;

	return
		(
			arg_dom_clust_id.segments
			? ( "[" + get_segments_string( arg_dom_clust_id.segments.get() ) + "]" )
			: "*"
		)
		+ (
			arg_include_cluster
			? ( "->" + to_string( arg_dom_clust_id.cluster_id ) )
			: ""
		);
}

/// \brief Generate a string describing the specified domain_cluster_ids
///
/// \relates domain_cluster_ids
string cath::clust::to_string(const domain_cluster_ids &arg_dom_clust_ids,  ///< The domain_cluster_ids to describe
                              const bool               &arg_include_cluster ///< Whether to include the cluster ID for each entry (default: true)
                              ) {
	return join(
		arg_dom_clust_ids
			| transformed( [&] (const domain_cluster_id &x) { return to_string( x, arg_include_cluster ); } ),
		", "
	);
}
