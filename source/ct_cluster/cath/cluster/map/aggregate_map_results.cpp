/// \map
/// \brief The aggregate_map_results class definitions

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

#include "aggregate_map_results.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/cluster/map/map_results.hpp"
#include "cath/cluster/new_cluster_data.hpp"
#include "cath/cluster/old_cluster_data.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/file/open_fstream.hpp"

using namespace cath::clust;
using namespace cath::common;

using boost::filesystem::path;
using boost::format;
using boost::lexical_cast;
using boost::numeric_cast;
using std::ofstream;
using std::string;

/// \brief Ctor from a clust_mapping_spec under which the mappings have been / will be performed
aggregate_map_results::aggregate_map_results(const clust_mapping_spec &prm_clust_mapping_spec ///< The clust_mapping_spec under which the mappings have been / will be performed
                                             ) noexcept : the_spec { prm_clust_mapping_spec } {
}

/// \brief Getter for whether this has had results added to it
const bool & aggregate_map_results::get_added_to() const {
	return added_to;
}

/// \brief Getter for the total number of old clusters encountered in the mapping
const size_t & aggregate_map_results::get_num_old_clusters() const {
	return num_old_clusters;
}

/// \brief Getter for the total number of new clusters encountered in the mapping
const size_t & aggregate_map_results::get_num_new_clusters() const {
	return num_new_clusters;
}

/// \brief Getter for the total number of clusters (on either side) that got mapped
const size_t & aggregate_map_results::get_num_mapped_clusters() const {
	return num_mapped_clusters;
}

/// \brief Getter for the total number of old entries encountered in the mapping
const size_t & aggregate_map_results::get_num_old_entries() const {
	return num_old_entries;
}

/// \brief Getter for the total number of new entries encountered in the mapping
const size_t & aggregate_map_results::get_num_new_entries() const {
	return num_new_entries;
}

/// \brief Getter for the total number of entries on either side that got mapped
const size_t & aggregate_map_results::get_num_mapped_entries() const {
	return num_mapped_entries;
}

/// \brief Getter for the number of entries which got a domain overlap of 0 *because there were no matching entries on the parent at all*
const size_t & aggregate_map_results::get_num_with_nothing_on_parent() const {
	return num_with_nothing_on_parent;
}

/// \brief Getter for the highest overlap fraction (over largest) for each of the old domains
const overlap_frac_distn & aggregate_map_results::get_highest_old_dom_overlap_fractions() const {
	return highest_old_dom_overlap_fractions;
}

/// \brief Getter for the highest overlap fraction for each of the old clusters
const overlap_frac_distn & aggregate_map_results::get_highest_old_clust_overlap_fractions() const {
	return highest_old_clust_overlap_fractions;
}


/// \brief Getter for the specification that was used to perform the mappings
const clust_mapping_spec & aggregate_map_results::get_clust_mapping_spec() const {
	return the_spec;
}

/// \brief Add the results from a new mapping
aggregate_map_results & aggregate_map_results::add_map_results(const map_results      &prm_map_results,  ///< The results to add
                                                               const old_cluster_data &prm_old_clusters, ///< The old clusters that were mapped from
                                                               const new_cluster_data &prm_new_clusters  ///< The new clusters that were mapped to
                                                               ) {
	if ( prm_map_results.the_spec != the_spec ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot add map_results that were performed under a different clust_mapping_spec to an aggregate_map_results"));
	}

	const size_t  curr_num_mapped_clusters = prm_map_results.chosen_maps.size();
	const size_t  curr_num_old_clusters    = get_num_clusters      ( prm_old_clusters );
	const size_t  curr_num_new_clusters    = get_num_clusters      ( prm_new_clusters );
	const size_t  curr_num_old_entries     = get_num_entries       ( prm_old_clusters );
	const size_t  curr_num_new_entries     = get_num_entries       ( prm_new_clusters );
	const size_t  curr_num_mapped_entries  = ::cath::clust::get_num_mapped_entries( prm_map_results  );

	if ( curr_num_mapped_clusters > curr_num_old_clusters || curr_num_mapped_clusters > curr_num_new_clusters ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Cannot have mapped more clusters than there are old/new clusters"));
	}

	num_old_clusters                    += curr_num_old_clusters;
	num_new_clusters                    += curr_num_new_clusters;
	num_mapped_clusters                 += curr_num_mapped_clusters;
	num_old_entries                     += curr_num_old_entries;
	num_new_entries                     += curr_num_new_entries;
	num_mapped_entries                  += curr_num_mapped_entries;
	num_with_nothing_on_parent          += prm_map_results.num_with_nothing_on_parent;

	highest_old_dom_overlap_fractions   += prm_map_results.highest_old_dom_overlap_fractions;
	highest_old_clust_overlap_fractions += prm_map_results.highest_old_clust_overlap_fractions;

	added_to                             = true;

	return *this;
}

/// \brief Make a aggregate_map_results containing one entry, the mapping between the specified clusters with the specified results
///
/// \relates aggregate_map_results
aggregate_map_results cath::clust::make_aggregate_map_results(const map_results      &prm_map_results,      ///< The results with which the aggregate_map_results should be initialised
                                                              const old_cluster_data &prm_old_cluster_data, ///< The old clusters that were mapped from
                                                              const new_cluster_data &prm_new_cluster_data  ///< The new clusters that were mapped to
                                                              ) {
	aggregate_map_results result{ prm_map_results.the_spec };
	result.add_map_results( prm_map_results, prm_old_cluster_data, prm_new_cluster_data );
	return result;
}

/// \brief Generate a Markdown summary of the specified map_results
///
/// \relates aggregate_map_results
string cath::clust::markdown_summary_string(const aggregate_map_results &prm_aggregate_map_results ///< The map_results to summarise
                                            ) {
	if ( ! prm_aggregate_map_results.get_added_to() ) {
		return "No mapping was performed";
	}

	const doub_vec percentile_list = { 25.0, 50.0, 75.0, 90.0, 95.0, 98.0, 99.0, 100.0 };

	const auto   &highest_old_dom_ol_fracs  = prm_aggregate_map_results.get_highest_old_dom_overlap_fractions();
	const auto   &highest_old_clst_ol_fracs = prm_aggregate_map_results.get_highest_old_clust_overlap_fractions();

	const size_t &num_old_clusters          = prm_aggregate_map_results.get_num_old_clusters();
	const size_t &num_new_clusters          = prm_aggregate_map_results.get_num_new_clusters();
	const size_t &num_mapped_clusters       = prm_aggregate_map_results.get_num_mapped_clusters();
	const size_t  num_unmapped_old_clusters = num_old_clusters - num_mapped_clusters;
	const size_t  num_unmapped_new_clusters = num_new_clusters - num_mapped_clusters;

	const size_t &num_old_entries           = prm_aggregate_map_results.get_num_old_entries();
	const size_t &num_new_entries           = prm_aggregate_map_results.get_num_new_entries();
	const size_t &num_mapped_entries        = prm_aggregate_map_results.get_num_mapped_entries();
	const size_t  num_unmapped_old_entries  = num_old_entries - num_mapped_entries;
	const size_t  num_unmapped_new_entries  = num_new_entries - num_mapped_entries;

	const double  pc_old_clusters_mapped    = 100.0 * numeric_cast<double>( num_mapped_clusters       ) / numeric_cast<double>( num_old_clusters );
	const double  pc_new_clusters_mapped    = 100.0 * numeric_cast<double>( num_mapped_clusters       ) / numeric_cast<double>( num_new_clusters );
	const double  pc_old_clusters_unmapped  = 100.0 * numeric_cast<double>( num_unmapped_old_clusters ) / numeric_cast<double>( num_old_clusters );
	const double  pc_new_clusters_unmapped  = 100.0 * numeric_cast<double>( num_unmapped_new_clusters ) / numeric_cast<double>( num_new_clusters );

	const double  pc_old_entries_mapped     = 100.0 * numeric_cast<double>( num_mapped_entries        ) / numeric_cast<double>( num_old_entries  );
	const double  pc_new_entries_mapped     = 100.0 * numeric_cast<double>( num_mapped_entries        ) / numeric_cast<double>( num_new_entries  );

	const double  pc_old_entries_unmapped   = 100.0 * numeric_cast<double>( num_unmapped_old_entries  ) / numeric_cast<double>( num_old_entries  );
	const double  pc_new_entries_unmapped   = 100.0 * numeric_cast<double>( num_unmapped_new_entries  ) / numeric_cast<double>( num_new_entries  );

	const double  min_equiv_dom_ol_pc       = 100.0 * prm_aggregate_map_results.get_clust_mapping_spec().get_min_equiv_dom_ol();
	const double  min_equiv_clust_ol_pc     = 100.0 * prm_aggregate_map_results.get_clust_mapping_spec().get_min_equiv_clust_ol();

	const string  dom_pc_str                = ( format( "%.1f" ) % min_equiv_dom_ol_pc ).str();

	return R"(Domain Mapping
==

This section describes how well the map-from domains could be mapped to new domains (and vice versa).
The quality of a mapping between a pair of domains is defined as the percentage overlap over the longer
domain (ie the percentage of the longer domain's residues shared with the other domain).
In this run, the cut-off for defining domain-equivalence was )" + dom_pc_str + R"(%.


Domains from Map-From Clusters
--

| Category                                                                          | Number | Percentage |
|-----------------------------------------------------------------------------------|--------|------------|
| All                                                                               |)" + ( format( "%7d" ) % num_old_entries          ).str() + R"( |     100.0% |
| &nbsp; ...of which:                                                               |        |            |
| &nbsp; &bull; Equivalence-mapped     (ie )" + dom_pc_str + R"( < overlap                         ) |)" + ( format( "%7d" ) % num_mapped_entries       ).str() + R"( | )" + ( ( format( "%9.1f" ) % pc_old_entries_mapped   ).str() ) + R"(% |
| &nbsp; &bull; Not equivalence-mapped (ie                         overlap ≤ )" + dom_pc_str + R"( ) |)" + ( format( "%7d" ) % num_unmapped_old_entries ).str() + R"( | )" + ( ( format( "%9.1f" ) % pc_old_entries_unmapped ).str() ) + R"(% |

Domains from New Clusters
--

| Category                                                       | Number | Percentage |
|----------------------------------------------------------------|--------|------------|
| All                                                            |)" + ( format( "%7d" ) % num_new_entries          ).str() + R"( |     100.0% |
| &nbsp; ...of which:                                            |        |            |
| &nbsp; &bull; Equivalence-mapped     (ie overlap > )" + dom_pc_str + R"(      ) |)" + ( format( "%7d" ) % num_mapped_entries       ).str() + R"( | )" + ( ( format( "%9.1f" ) % pc_new_entries_mapped   ).str() ) + R"(% |
| &nbsp; &bull; Not equivalence-mapped (ie overlap ≤ )" + dom_pc_str + R"(      ) |)" + ( format( "%7d" ) % num_unmapped_new_entries ).str() + R"( | )" + ( ( format( "%9.1f" ) % pc_new_entries_unmapped ).str() ) + R"(% |

Distribution of Domain Mapping Percentages for Domains from Map-From Clusters
--

Including completely-unmapped domains:

)"
	+ percentile_markdown_table(
		highest_old_dom_ol_fracs,
		percentile_list,
		"Percentile through distribution of mapping percentages",
		"Mapping Percentage",
		overlap_frac_distn::zeroes_policy::INCLUDE
	) + R"(

Excluding completely-unmapped domains:

)"
	+ percentile_markdown_table(
		highest_old_dom_ol_fracs,
		percentile_list,
		"Percentile through distribution of mapping percentages",
		"Mapping Percentage",
		overlap_frac_distn::zeroes_policy::EXCLUDE
	) + R"(

Histogram of domain mapping percentages

)"
	+ histogram_markdown_table(
		highest_old_dom_ol_fracs,
		"Domain overlap",
		"Number of map-from domains",
		"Percentage of map-from domains",
		prm_aggregate_map_results.get_num_with_nothing_on_parent()
	)
	+ R"(


Cluster Mapping
==

This section describes how well the old clusters could be mapped to the new clusters (and vice versa).
The quality of a mapping between a pair of clusters is defined as the percentage of the domains in
the map-from cluster that have an equivalent domain in the new cluster.
In this run, the cutoff for defining cluster-equivalence was )" + ( format( "%.1f" ) % min_equiv_clust_ol_pc ).str() + R"(%.

Also, for clusters to be considered equivalents, the map-from cluster's members must
be equivalent to > )"
	+ lexical_cast<string>( 100.0 * clust_mapping_spec::MIN_EQUIV_FRAC_OF_NEW_CLUST        )
	+ R"(% of the working cluster's entries and > )"
	+ lexical_cast<string>( 100.0 * clust_mapping_spec::MIN_EQUIV_FRAC_OF_NEW_CLUST_EQUIVS )
	+ R"(% of those that have an equivalence.

Map-From Clusters
--

| Category                             | Number | Percentage |
|--------------------------------------|--------|------------|
| All                                  |)" + ( format( "%7d" ) % num_old_clusters          ).str() + R"( |     100.0% |
| &nbsp; ...of which:                  |        |            |
| &nbsp; &bull; Equivalence-mapped     |)" + ( format( "%7d" ) % num_mapped_clusters       ).str() + R"( | )" + ( ( format( "%9.1f" ) % pc_old_clusters_mapped   ).str() ) + R"(% |
| &nbsp; &bull; Not equivalence-mapped |)" + ( format( "%7d" ) % num_unmapped_old_clusters ).str() + R"( | )" + ( ( format( "%9.1f" ) % pc_old_clusters_unmapped ).str() ) + R"(% |

New Clusters
--

| Category                             | Number | Percentage |
|--------------------------------------|--------|------------|
| All                                  |)" + ( format( "%7d" ) % num_new_clusters          ).str() + R"( |     100.0% |
| &nbsp; ...of which:                  |        |            |
| &nbsp; &bull; Equivalence-mapped     |)" + ( format( "%7d" ) % num_mapped_clusters       ).str() + R"( | )" + ( ( format( "%9.1f" ) % pc_new_clusters_mapped   ).str() ) + R"(% |
| &nbsp; &bull; Not equivalence-mapped |)" + ( format( "%7d" ) % num_unmapped_new_clusters ).str() + R"( | )" + ( ( format( "%9.1f" ) % pc_new_clusters_unmapped ).str() ) + R"(% |

Distribution of Cluster Mapping Percentages for Map-From Clusters
--

Including completely-unmapped clusters:

)"
	+ percentile_markdown_table(
		highest_old_clst_ol_fracs,
		percentile_list,
		"Percentile through distribution of mapping percentages",
		"Mapping Percentage",
		overlap_frac_distn::zeroes_policy::INCLUDE
	)
	+ R"(

Excluding completely-unmapped clusters:

)"
	+ percentile_markdown_table(
		highest_old_clst_ol_fracs,
		percentile_list,
		"Percentile through distribution of mapping percentages",
		"Mapping Percentage",
		overlap_frac_distn::zeroes_policy::EXCLUDE
	)
	+ "\n\n";
}

/// TODO: Consider splitting "Domains from Map-From Clusters" numbers entry:
///  * Not equivalence-mapped (ie     overlap ≤ X )
/// into:
///  * Insufficiently-mapped  (ie 0 < overlap ≤ X )
///  * Completely-unmapped    (ie     overlap = 0 )

/// \brief Write a Markdown summary of the specified aggregate_map_results to the specified file
///
/// \relates aggregate_map_results
void cath::clust::write_markdown_summary_string_to_file(const path                  &prm_output_file,          ///< The file to which the Markdown summary should be written
                                                        const aggregate_map_results &prm_aggregate_map_results ///< The aggregate_map_results to summarise
                                                        ) {
	ofstream md_ostream;
	open_ofstream( md_ostream, prm_output_file );
	md_ostream << markdown_summary_string( prm_aggregate_map_results );
	md_ostream.close();
}
