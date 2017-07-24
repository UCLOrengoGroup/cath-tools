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

#include "common/file/slurp.hpp"
#include "test/global_test_constants.hpp"

using namespace cath::common;
using namespace cath::test;

using boost::filesystem::path;
using std::string;

/// \brief Test constant for the cath-resolve-hits test data directory
path map_clusters_fixture::map_clusters_test_data_dir() {
	return global_test_constants::TEST_SOURCE_DATA_DIR() / "map-clusters" ;
}

/// \brief Get the test file for the cath-map-clusters example input file
path map_clusters_fixture::eg_input_file() {
	return map_clusters_test_data_dir() / "3.30.910.10.later";
}

/// \brief Get the test file for the cath-map-clusters example input map-from file
path map_clusters_fixture::eg_input_mapfrom_file() {
	return map_clusters_test_data_dir() / "3.30.910.10.v4_0_0";
}

/// \brief Get the test file for the cath-map-clusters example map-from result file
path map_clusters_fixture::eg_mapfrom_result_file() {
	return map_clusters_test_data_dir() / "3.30.910.10.result.mapfrom";
}

/// \brief Get the test file for the cath-map-clusters example renumber only_result file
path map_clusters_fixture::eg_renumber_only_result_file() {
	return map_clusters_test_data_dir() / "3.30.910.10.result.renumber_only";
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
