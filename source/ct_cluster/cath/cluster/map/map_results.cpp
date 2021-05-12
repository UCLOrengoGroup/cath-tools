/// \map
/// \brief The map_results class definitions

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

#include "map_results.hpp"

#include <filesystem>
#include <fstream>

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/numeric.hpp>

#include "cath/cluster/map/aggregate_map_results.hpp"
#include "cath/cluster/new_cluster_data.hpp"
#include "cath/cluster/old_cluster_data.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/file/open_fstream.hpp"

using namespace ::cath::common;
using namespace ::std::literals::string_literals;

using ::boost::adaptors::transformed;
using ::boost::algorithm::join;
using ::boost::format;
using ::std::filesystem::path;
using ::std::ofstream;
using ::std::string;

//               new_0[25]  new_1[23]  new_2[ 2]  new_3[64]  new_4[37]  new_5[33]
//  old_0 [20] |          |          |          |     18   |          |          |
//  old_1 [23] |          |          |          |          |          |     17   |
//  old_2 [15] |     14   |          |          |          |          |          |
//  old_3 [22] |          |     16   |          |          |          |          |
//  old_4 [ 2] |          |          |      2   |          |          |          |
//  old_5 [52] |          |          |          |      2   |     28   |          |
//
//  Mapped  18 between old_0:  20 [ 18 mapped] and new_3:  64 [ 20 mapped]   90.0%  100.0%   28.1%   90.0%
//  Mapped  17 between old_1:  23 [ 17 mapped] and new_5:  33 [ 17 mapped]   73.9%  100.0%   51.5%  100.0%
//  Mapped  14 between old_2:  15 [ 14 mapped] and new_0:  25 [ 14 mapped]   93.3%  100.0%   56.0%  100.0%
//  Mapped  16 between old_3:  22 [ 16 mapped] and new_1:  23 [ 16 mapped]   72.7%  100.0%   69.6%  100.0%
//  Mapped   2 between old_4:   2 [  2 mapped] and new_2:   2 [  2 mapped]  100.0%  100.0%  100.0%  100.0%
//  Mapped   2 between old_5:  52 [ 30 mapped] and new_3:  64 [ 20 mapped]    3.8%    6.7%    3.1%   10.0%
//  Mapped  28 between old_5:  52 [ 30 mapped] and new_4:  37 [ 28 mapped]   53.8%   93.3%   75.7%  100.0%

/// \brief Get the name of the new, unmapped cluster of the specified index
///        based on the specified last preceding index (or 0 if not mapping from
///        or none if mapping from clusters that don't all have numeric names)
string cath::clust::detail::get_name_of_new_unmapped_cluster_of_index(const boost::optional<ptrdiff_t> &prm_last_preceding_index, ///< The last preceding index (or 0 if not mapping from or none if mapping from clusters that don't all have numeric names)
                                                                      const size_t                     &prm_index                 ///< The index of the new, unmapped cluster
                                                                      ) {
	using ::std::to_string;
	if ( prm_last_preceding_index ) {
		return to_string( *prm_last_preceding_index + 1l + debug_numeric_cast<ptrdiff_t>( prm_index ) );
	}
	return "new_cmc_cluster_" + to_string( prm_index + 1 ) ;
}

/// \brief Get the number of mapped entries implied by the specified map_results
///
/// \pre The total mapped from the old clusters must match the total mapped to the new clusters
///      else an invalid_argument_exception will be thrown
///
/// \relates map_results
size_t cath::clust::get_num_mapped_entries(const map_results &prm_map_results ///< The map_results to query
                                           ) {
	const size_t num_mapped_olds = boost::accumulate( prm_map_results.num_mapped_by_old_cluster, 0_z );
	const size_t num_mapped_news = boost::accumulate( prm_map_results.num_mapped_by_new_cluster, 0_z );

	if ( num_mapped_olds != num_mapped_news ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot retrieve the number of mapped entries if the answer is inconsistent for old/new"));
	}
	return num_mapped_olds;
}

/// \brief Generate a string describing the specified map_results
///
/// \relates map_results
string cath::clust::results_string(const old_cluster_data_opt &prm_old_clusters,   ///< The old clusters
                                   const new_cluster_data     &prm_new_clusters,   ///< The new clusters
                                   const map_results          &prm_map_results,    ///< The map_results
                                   const str_opt              &prm_batch_id,       ///< An optional identifier for the mapped batch
                                   const bool                 &prm_include_headers ///< Whether to include the headers for the columns
                                   ) {
	const auto &chosen_maps                  = prm_map_results.chosen_maps;
	const auto &other_maps                   = prm_map_results.other_maps;
	const auto &unmapped_new_cluster_indices = prm_map_results.unmapped_new_cluster_indices;

	// If not mapping, start from 1
	// If mapping and numeric, start from last + 1
	// If mapping and not numeric
	const auto precede_index = largest_number_if_names_all_numeric_integers_of_val_if_none( prm_old_clusters, 0 );

	if ( ! prm_old_clusters && ( ! chosen_maps.empty() || ! other_maps.empty() ) ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Argh"));
	}

	return (
			prm_include_headers
				? "# cluster-id suggested-name"
					+ ( prm_batch_id ? " batch-id"s : ""s )
					+ "\n"
				: ""
		)
		+ join(
			chosen_maps
				| transformed( [&] (const potential_map &the_map) {
					return
						  get_name_of_cluster_of_id(  prm_new_clusters, the_map.new_cluster_idx )
						+ " "
						+ get_name_of_cluster_of_id( *prm_old_clusters, the_map.old_cluster_idx )
						+ ( prm_batch_id ? ( " " + *prm_batch_id ) : "" )
						+ "\n";
				} ),
			""
		)
		// + join(
		// 	other_maps
		// 		| transformed( [&] (const potential_map &the_map) {
		// 			return
		// 				  get_name_of_cluster_of_id(  prm_new_clusters, the_map.new_cluster_idx )
		// 				+ " "
		// 				+ get_name_of_cluster_of_id( *prm_old_clusters, the_map.old_cluster_idx )
		// 				+ ( prm_batch_id ? ( " " + *prm_batch_id ) : "" )
		// 				+ " other\n";
		// 		} ),
		// 	""
		// )
		+ join(
			indices( unmapped_new_cluster_indices.size() )
				| transformed( [&] (const size_t &new_cluster_index_index) {
					const size_t &new_cluster_index = unmapped_new_cluster_indices[ new_cluster_index_index ];
					return
						  get_name_of_cluster_of_id( prm_new_clusters,  new_cluster_index )
						+ " "
						+ detail::get_name_of_new_unmapped_cluster_of_index( precede_index, new_cluster_index_index )
						+ ( prm_batch_id ? ( " " + *prm_batch_id ) : "" )
						+ "\n";
				} ),
			""
		);
}

/// \brief Generate a string describing the specified map_results
///
/// \relates map_results
string cath::clust::longer_results_string(const old_cluster_data_opt &prm_old_clusters, ///< The old clusters
                                          const new_cluster_data     &prm_new_clusters, ///< The new clusters
                                          const map_results          &prm_map_results,  ///< The map_results
                                          const str_opt              &prm_batch_id      ///< An optional identifier for the mapped batch
                                          ) {
	const auto &chosen_maps                  = prm_map_results.chosen_maps;
	const auto &other_maps                   = prm_map_results.other_maps;
	const auto &num_mapped_by_new_cluster    = prm_map_results.num_mapped_by_new_cluster;
	const auto &num_mapped_by_old_cluster    = prm_map_results.num_mapped_by_old_cluster;
	const auto &unmapped_new_cluster_indices = prm_map_results.unmapped_new_cluster_indices;

	// If not mapping, start from 1
	// If mapping and numeric, start from last + 1
	// If mapping and not numeric
	const auto precede_index = largest_number_if_names_all_numeric_integers_of_val_if_none( prm_old_clusters, 0 );

	if ( ! prm_old_clusters && ( ! chosen_maps.empty() || ! other_maps.empty() ) ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Argh"));
	}

	const auto map_summary_fn = [&] (const potential_map &the_map) {
		const size_t &old_cluster_idx   = the_map.old_cluster_idx;
		const size_t &new_cluster_idx   = the_map.new_cluster_idx;
		const size_t &num_mapped        = the_map.num_mapped;
		const size_t &num_in_new        = get_size_of_cluster_of_id( prm_new_clusters, new_cluster_idx );
		const size_t &num_in_old        = num_entries( ( *prm_old_clusters ) [ old_cluster_idx ] );
		const size_t &num_mapped_in_new = num_mapped_by_new_cluster[ new_cluster_idx ];
		const size_t &num_mapped_in_old = num_mapped_by_old_cluster[ old_cluster_idx ];

		const string new_name = ( format( "%3d" ) % get_name_of_cluster_of_id(  prm_new_clusters, new_cluster_idx ) ).str();
		const string old_name = ( format( "%3d" ) % get_name_of_cluster_of_id( *prm_old_clusters, old_cluster_idx ) ).str();

		return
			  new_name
			+ " to "
			+ old_name
			+ " [share "
			+ ( format( "%3d" ) % num_mapped        ).str()
			+ " equivs ie: "
			+ ( format( "%5.1f" ) % ( static_cast<double>( num_mapped ) * 100.0 / static_cast<double>( num_in_new        ) ) ).str()
			+ R"(% of )"
			+ new_name
			+ "'s "
			+ ( format( "%3d" ) % num_in_new        ).str()
			+ " members and  "
			+ ( format( "%5.1f" ) % ( static_cast<double>( num_mapped ) * 100.0 / static_cast<double>( num_mapped_in_new ) ) ).str()
			+ R"(% of its )"
			+ ( format( "%3d" ) % num_mapped_in_new ).str()
			+ " equivs; "
			+ ( format( "%5.1f" ) % ( static_cast<double>( num_mapped ) * 100.0 / static_cast<double>( num_in_old        ) ) ).str()
			+ R"(% of )"
			+ old_name
			+ "'s "
			+ ( format( "%3d" ) % num_in_old        ).str()
			+ "  members and "
			+ ( format( "%5.1f" ) % ( static_cast<double>( num_mapped ) * 100.0 / static_cast<double>( num_mapped_in_old ) ) ).str()
			+ R"(% of its )"
			+ ( format( "%3d" ) % num_mapped_in_old ).str()
			+ " equivs]";
	};

	return
		join(
			chosen_maps
				| transformed( [&] (const potential_map &the_map) {
					return "RENAME      " + ( prm_batch_id ? *prm_batch_id : "" ) + map_summary_fn( the_map );
				} ),
			"\n"
		)
		+ "\n"
		+ join(
			indices( unmapped_new_cluster_indices.size() )
				| transformed( [&] (const size_t &new_cluster_index_index) {
					const size_t &new_cluster_index = unmapped_new_cluster_indices[ new_cluster_index_index ];
					return
						  "RENAME      "
						+ ( prm_batch_id ? " " +*prm_batch_id : "" )
						+ get_name_of_cluster_of_id( prm_new_clusters,  new_cluster_index )
						+ " to "
						+ detail::get_name_of_new_unmapped_cluster_of_index( precede_index, new_cluster_index_index );
				} ),
			"\n"
		)
		+ "\n"
		+ join(
			other_maps
				| transformed( [&] (const potential_map &the_map) {
					return "#ALSOLINKED " + map_summary_fn( the_map );
				} ),
			"\n"
		)
		+ "\n";
}

/// \brief Generate a Markdown summary of the specified map_results
///
/// \relates map_results
string cath::clust::markdown_summary_string(const old_cluster_data_opt &prm_old_clusters, ///< The old clusters
                                            const new_cluster_data     &prm_new_clusters, ///< The new clusters
                                            const map_results          &prm_map_results   ///< The map_results to summarise
                                             ) {
	using ::std::to_string;

	if ( ! prm_old_clusters ) {
		return "Mapping was not performed\n";
	}

	return markdown_summary_string( make_aggregate_map_results(
		prm_map_results,
		*prm_old_clusters,
		prm_new_clusters
	) );
}

/// \brief Write a Markdown summary of the specified map_results to the specified file
///
/// \relates map_results
void cath::clust::write_markdown_summary_string_to_file(const path                 &prm_output_file,  ///< The file to which the Markdown summary should be written
                                                        const old_cluster_data_opt &prm_old_clusters, ///< The old clusters
                                                        const new_cluster_data     &prm_new_clusters, ///< The new clusters
                                                        const map_results          &prm_map_results   ///< The map_results to summarise
                                                        ) {
	ofstream md_ostream = open_ofstream( prm_output_file );
	md_ostream << markdown_summary_string( prm_old_clusters, prm_new_clusters, prm_map_results );
	md_ostream.close();
}
