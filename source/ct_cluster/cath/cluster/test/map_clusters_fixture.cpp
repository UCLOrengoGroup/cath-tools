/// \file
/// \brief The map_clusters_fixture class definitions

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

#include "map_clusters_fixture.hpp"

#include "cath/common/file/slurp.hpp"
#include "cath/test/global_test_constants.hpp"

using namespace cath::common;
using namespace cath::test;

using boost::filesystem::path;
using std::string;

/// \brief Test constant for the cath-map-clusters test data directory
path map_clusters_fixture::map_clusters_test_data_dir() {
	return global_test_constants::TEST_SOURCE_DATA_DIR() / "map-clusters" ;
}

/// \brief Test constant for the cath-map-clusters to_clust_ol_thresholds test data directory
path map_clusters_fixture::to_clust_ol_thresholds_test_data_dir() {
	return map_clusters_test_data_dir() / "to_clust_ol_thresholds" ;
}

/// \brief Test constant for the cath-map-clusters quirky_cluster_memberships test data directory
path map_clusters_fixture::eg_quirky_clustmemb_test_data_dir() {
	return map_clusters_test_data_dir() / "quirky_cluster_memberships" ;
}

/// \brief Test constant for the cath-map-clusters to_clust_ol_thresholds input file
path map_clusters_fixture::to_clust_ol_thresholds_input_file() {
	return to_clust_ol_thresholds_test_data_dir() / "to_clust_ol_thresholds.mapto" ;
}

/// \brief Test constant for the cath-map-clusters to_clust_ol_thresholds map-from file
path map_clusters_fixture::to_clust_ol_thresholds_mapfrom_file() {
	return to_clust_ol_thresholds_test_data_dir() / "to_clust_ol_thresholds.mapfrom" ;
}

/// \brief Test constant for the cath-map-clusters to_clust_ol_thresholds result file
path map_clusters_fixture::to_clust_ol_thresholds_result_file() {
	return to_clust_ol_thresholds_test_data_dir() / "to_clust_ol_thresholds.result" ;
}

/// \brief Get the test file for the cath-map-clusters help usage file
path map_clusters_fixture::help_usage_file() {
	return map_clusters_test_data_dir() / "help_usage";
}

/// \brief Get the test file for the cath-map-clusters empty result file
path map_clusters_fixture::empty_result_file() {
	return map_clusters_test_data_dir() / "result.empty";
}

/// \brief Get the test file for the cath-map-clusters example input file
path map_clusters_fixture::eg_input_file() {
	return map_clusters_test_data_dir() / "example.later";
}

/// \brief Get the test file for the cath-map-clusters example input map-from file
path map_clusters_fixture::eg_input_mapfrom_file() {
	return map_clusters_test_data_dir() / "example.v4_0_0";
}

/// \brief Get the test file for the cath-map-clusters example map-from result file
path map_clusters_fixture::eg_mapfrom_result_file() {
	return map_clusters_test_data_dir() / "example.result.mapfrom";
}

/// \brief Get the test file for the cath-map-clusters example renumber only result file
path map_clusters_fixture::eg_renumber_only_result_file() {
	return map_clusters_test_data_dir() / "example.result.renumber_only";
}

/// \brief Get the test file for the cath-map-clusters example append batch id result file
path map_clusters_fixture::eg_append_batch_id_result_file() {
	return map_clusters_test_data_dir() / "example.result.append_batch_id";
}

/// \brief Get the test file for the cath-map-clusters example mapping summary file
path map_clusters_fixture::eg_summary_mapping_file() {
	return map_clusters_test_data_dir() / "example.summary.mapping";
}

/// \brief Get the test file for the cath-map-clusters example renumbering summary file
path map_clusters_fixture::eg_renumbering_summary_file() {
	return map_clusters_test_data_dir() / "example.summary.renumbering";
}

/// \brief Get the test file for the cath-map-clusters example batch input file
path map_clusters_fixture::eg_batch_input_file() {
	return map_clusters_test_data_dir() / "batch_file";
}

/// \brief Get the test file for the cath-map-clusters example batch result file
path map_clusters_fixture::eg_batch_result_file() {
	return map_clusters_test_data_dir() / "batch_result";
}

/// \brief Get the test file for the cath-map-clusters example entry-mapping result file
path map_clusters_fixture::eg_entry_mapping_result_file() {
	return map_clusters_test_data_dir() / "entry_mapping_result";
}

/// \brief Get the test file for the cath-map-clusters example simple map-from result file
path map_clusters_fixture::eg_simple_mapfrom_result_file() {
	return map_clusters_test_data_dir() / "example.result.simplemapfrom";
}

/// \brief Get the test file for the cath-map-clusters example map-from result file with dom_ol_50
path map_clusters_fixture::eg_mapfrom_dom_ol_50_result_file() {
	return map_clusters_test_data_dir() / "example.result.mapfrom.dom_ol_50";
}

/// \brief Get the test file for the cath-map-clusters example map-from result file with clust_ol_50
path map_clusters_fixture::eg_mapfrom_clust_ol_50_result_file() {
	return map_clusters_test_data_dir() / "example.result.mapfrom.clust_ol_50";
}

/// \brief Get the test file for the cath-map-clusters example map-from result file with dom_ol_100
path map_clusters_fixture::eg_mapfrom_dom_ol_100_result_file() {
	return map_clusters_test_data_dir() / "example.result.mapfrom.dom_ol_100";
}

/// \brief Get the test file for the cath-map-clusters example map-from result file with clust_ol_100
path map_clusters_fixture::eg_mapfrom_clust_ol_100_result_file() {
	return map_clusters_test_data_dir() / "example.result.mapfrom.clust_ol_100";
}

/// \brief Get the test file for the cath-map-clusters example input file with non-numeric cluster names
path map_clusters_fixture::eg_input_non_numeric_file() {
	return map_clusters_test_data_dir() / "example.cluster_membership.non_numeric_names";
}

/// \brief Get the test file for the cath-map-clusters example result file after using eg_input_non_numeric_file() as the from file
path map_clusters_fixture::eg_input_non_numeric_fromresult() {
	return map_clusters_test_data_dir() / "example.cluster_membership.non_numeric_names.fromresult";
}

/// \brief Get the test file for the cath-map-clusters example result file after using eg_input_non_numeric_file() as the to file
path map_clusters_fixture::eg_input_non_numeric_toresult() {
	return map_clusters_test_data_dir() / "example.cluster_membership.non_numeric_names.toresult";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting clashing segments
path map_clusters_fixture::eg_input_clashing_segments_file() {
	return eg_quirky_clustmemb_test_data_dir() / "clashing_segments";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting clashing segments with different names
path map_clusters_fixture::eg_input_clashing_segments_w_diff_names_file() {
	return eg_quirky_clustmemb_test_data_dir() / "clashing_segments_w_diff_names";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting mixed WCDs and segments
path map_clusters_fixture::eg_input_mixed_wcds_and_segments_file() {
	return eg_quirky_clustmemb_test_data_dir() / "mixed_wcds_and_segments";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting mixed WCDs and segments within a cluster
path map_clusters_fixture::eg_input_mixed_wcds_and_segments_within_cluster_file() {
	return eg_quirky_clustmemb_test_data_dir() / "mixed_wcds_and_segments_within_cluster";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting repeated segments
path map_clusters_fixture::eg_input_repeated_segments_file() {
	return eg_quirky_clustmemb_test_data_dir() / "repeated_segments";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting repeated segments with different names
path map_clusters_fixture::eg_input_repeated_segments_w_diff_names_file() {
	return eg_quirky_clustmemb_test_data_dir() / "repeated_segments_w_diff_names";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting repeated WCDs
path map_clusters_fixture::eg_input_repeated_wcds_file() {
	return eg_quirky_clustmemb_test_data_dir() / "repeated_wcds";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting backward segment
path map_clusters_fixture::eg_input_backward_segment() {
	return eg_quirky_clustmemb_test_data_dir() / "backward_segment";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting misordered segments
path map_clusters_fixture::eg_input_misordered_segments() {
	return eg_quirky_clustmemb_test_data_dir() / "misordered_segments";
}

/// \brief Get the test file for the cath-map-clusters example input file exhibiting zero start
path map_clusters_fixture::eg_input_zero_start() {
	return eg_quirky_clustmemb_test_data_dir() / "zero_start";
}



/// \brief Get the test string for the cath-map-clusters example input file
string map_clusters_fixture::eg_input_str() {
	return slurp( eg_input_file() );
}

/// \brief Get the test string for the cath-map-clusters example input map-from file
string map_clusters_fixture::eg_input_mapfrom_str() {
	return slurp( eg_input_mapfrom_file() );
}

/// \brief Get the test string for the cath-map-clusters example map-from result file
string map_clusters_fixture::eg_mapfrom_result_str() {
	return slurp( eg_mapfrom_result_file() );
}

/// \brief Get the test string for the cath-map-clusters example renumber only_result file
string map_clusters_fixture::eg_renumber_only_result_str() {
	return slurp( eg_renumber_only_result_file() );
}
