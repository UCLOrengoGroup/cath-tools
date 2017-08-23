/// \file
/// \brief The new_cluster_data class definitions

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

#include "new_cluster_data.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "common/algorithm/sort_uniq_build.hpp"

using namespace cath::common;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::irange;
using std::pair;
using std::string;

/// \brief Generate a string describing the specified new_cluster_data
///
/// \relates new_cluster_data
std::string cath::clust::to_string(const new_cluster_data &arg_new_cluster_data ///< The new_cluster_data to describe
                                   ) {
	using std::to_string;

	// Grab the id_of_seq_name of the new_cluster_data
	const id_of_str_bidirnl &id_of_seq_name = get_id_of_seq_name( arg_new_cluster_data );

	// Get a sorted list of the sequence names
	const auto sorted_seq_names = sort_build<str_vec>( id_of_seq_name );

	// Return a string summarising the new_cluster_data
	return "new_cluster_data["
		+ to_string( get_num_clusters( arg_new_cluster_data ) )
		+ " clusters, cluster_sizes{ "
		+ join(
			irange( 0_z, get_num_clusters( arg_new_cluster_data ) )
				| transformed( [&] (const size_t &x) {
					return to_string( x )
						+ R"((")"
						+ get_name_of_cluster_of_id( arg_new_cluster_data, x )
						+ R"("):)"
						+ to_string( get_size_of_cluster_of_id( arg_new_cluster_data, x ) );
				} ),
			", "
		)
		+ " }, cluster_index_by_seq_regions{ "
		+ join(
			sorted_seq_names
				| filtered( [&] (const string &x) {
					return has_domain_cluster_ids_of_seq_name( arg_new_cluster_data, x );
				} )
				| transformed( [&] (const string &x) {
					return x + ":(" + to_string( get_domain_cluster_ids_of_seq_name( arg_new_cluster_data, x ) ) + ")";
				} ),
			", "
		)
		+ " } ]";
}
