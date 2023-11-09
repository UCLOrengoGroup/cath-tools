/// \file
/// \brief The map_clusters_fixture header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_TEST_MAP_CLUSTERS_FIXTURE_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_TEST_MAP_CLUSTERS_FIXTURE_HPP

#include <filesystem>
#include <string>

namespace cath::test {

	/// \brief Store data that can be reused for multiple map_clusters tests
	class map_clusters_fixture {
	protected:
		static ::std::filesystem::path map_clusters_test_data_dir();
		static ::std::filesystem::path to_clust_ol_thresholds_test_data_dir();
		static ::std::filesystem::path eg_quirky_clustmemb_test_data_dir();

		static ::std::filesystem::path to_clust_ol_thresholds_input_file();
		static ::std::filesystem::path to_clust_ol_thresholds_mapfrom_file();
		static ::std::filesystem::path to_clust_ol_thresholds_result_file();

		static ::std::filesystem::path help_usage_file();

		static ::std::filesystem::path empty_result_file();

		static ::std::filesystem::path eg_input_file();
		static ::std::filesystem::path eg_input_mapfrom_file();
		static ::std::filesystem::path eg_mapfrom_result_file();
		static ::std::filesystem::path eg_renumber_only_result_file();
		static ::std::filesystem::path eg_append_batch_id_result_file();
		static ::std::filesystem::path eg_summary_mapping_file();
		static ::std::filesystem::path eg_renumbering_summary_file();
		static ::std::filesystem::path eg_batch_input_file();
		static ::std::filesystem::path eg_batch_result_file();
		static ::std::filesystem::path eg_entry_mapping_result_file();

		static ::std::filesystem::path eg_simple_mapfrom_result_file();
		static ::std::filesystem::path eg_mapfrom_dom_ol_50_result_file();
		static ::std::filesystem::path eg_mapfrom_clust_ol_50_result_file();
		static ::std::filesystem::path eg_mapfrom_dom_ol_100_result_file();
		static ::std::filesystem::path eg_mapfrom_clust_ol_100_result_file();

		static ::std::filesystem::path eg_input_non_numeric_file();
		static ::std::filesystem::path eg_input_non_numeric_fromresult();
		static ::std::filesystem::path eg_input_non_numeric_toresult();

		static ::std::filesystem::path eg_input_clashing_segments_file();
		static ::std::filesystem::path eg_input_clashing_segments_w_diff_names_file();
		static ::std::filesystem::path eg_input_mixed_wcds_and_segments_file();
		static ::std::filesystem::path eg_input_mixed_wcds_and_segments_within_cluster_file();
		static ::std::filesystem::path eg_input_repeated_segments_file();
		static ::std::filesystem::path eg_input_repeated_segments_w_diff_names_file();
		static ::std::filesystem::path eg_input_repeated_wcds_file();
		static ::std::filesystem::path eg_input_backward_segment();
		static ::std::filesystem::path eg_input_misordered_segments();
		static ::std::filesystem::path eg_input_zero_start();

		static std::string eg_input_str();
		static std::string eg_input_mapfrom_str();
		static std::string eg_mapfrom_result_str();
		static std::string eg_renumber_only_result_str();

	};

} // namespace cath::test

#endif // CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_TEST_MAP_CLUSTERS_FIXTURE_HPP
