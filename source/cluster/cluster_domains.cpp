/// \file
/// \brief The cluster_domains class definitions

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

#include "cluster_domains.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "common/container/id_of_str_bidirnl.hpp"

using namespace cath::clust::detail;
using namespace cath::common;

using boost::accumulate;
using boost::adaptors::transformed;
using boost::algorithm::join;
using std::string;

/// \brief Get the number of entries in the specified cluster_domains
///
/// \relates cluster_domains
size_t cath::clust::num_entries(const cluster_domains &arg_cluster_domains ///< The cluster_domains to query
                                ) {
	return accumulate(
		arg_cluster_domains
			| transformed( [] (const seq_id_and_domain_cluster_ids_pair &x) {
				return x.dom_cluster_ids.size();
			} ),
		0_z
	);
}

/// \brief Generate a string describing the specified cluster_domains
///
/// \relates cluster_domains
string cath::clust::to_string(const cluster_domains   &arg_cluster_domains, ///< The cluster_domains to describe
                              const id_of_str_bidirnl &arg_ider             ///< The id_of_str_bidirnl that mapped from seq name to ID
                              ) {
	return join(
		arg_cluster_domains
			| transformed( [&] (const seq_id_and_domain_cluster_ids_pair &x) {
				return
					  arg_ider.get_name_of_id( x.seq_id )
					+ "("
					+ to_string( x.dom_cluster_ids, false )
					+ ")";
			} ),
		", "
	);
}
