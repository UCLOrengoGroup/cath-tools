/// \file
/// \brief The old_cluster_data class definitions

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

#include "old_cluster_data.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "common/boost_addenda/range/indices.hpp"

using namespace cath::clust;
using namespace cath::common;

using boost::adaptors::transformed;
using boost::optional;
using boost::algorithm::join;

/// \brief Generate a string describing the specified old_cluster_data
///
/// \relates old_cluster_data
std::string cath::clust::to_string(const old_cluster_data &arg_old_cluster_data ///< The old_cluster_data to describe
                                   ) {
	using std::to_string;

	// Return a string summarising the old_cluster_data
	return "old_cluster_data["
		+ to_string( get_num_clusters( arg_old_cluster_data ) )
		+ " clusters, clusters{ "
		+ join(
			indices( get_num_clusters( arg_old_cluster_data ) )
				| transformed( [&] (const size_t &x) {
					return to_string( x )
						+ R"((")"
						+ get_name_of_cluster_of_id( arg_old_cluster_data, x )
						+ R"("): )"
						+ to_string(
							arg_old_cluster_data[ x ],
							arg_old_cluster_data.get_id_of_seq_name()
						);
				} ),
			", "
		)
		+ " } ]";
}

/// \brief Get the largest number if the names are all numeric or none otherwise
///
/// \relates old_cluster_data
optional<ptrdiff_t> cath::clust::largest_number_if_names_all_numeric_integers(const old_cluster_data &arg_old_cluster_data ///< The old_cluster_data to describe
                                                                              ) {
	return largest_number_if_names_all_numeric_integers(
		arg_old_cluster_data.get_clust_info().get_ider()
	);
}

/// \brief Get the largest number if the names are all numeric
///        or arg_value if no old_cluster_data is present
///        or none if the cluster names aren't all numeric
///
/// \relates old_cluster_data
optional<ptrdiff_t> cath::clust::largest_number_if_names_all_numeric_integers_of_val_if_none(const old_cluster_data_opt &arg_old_cluster_data, ///< The old_cluster_data to describe
                                                                                             const ptrdiff_t            &arg_value             ///< The value to return if no old_cluster_data is present
                                                                                             ) {
	return arg_old_cluster_data ? largest_number_if_names_all_numeric_integers( *arg_old_cluster_data )
	                            : arg_value;
}