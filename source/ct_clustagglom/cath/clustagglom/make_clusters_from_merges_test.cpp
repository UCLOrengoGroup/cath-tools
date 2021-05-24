/// \file
/// \brief The make_clusters_from_merges test suite

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

#include <filesystem>

#include <boost/test/unit_test.hpp>

#include "cath/clustagglom/clustagglom_fixture.hpp"
#include "cath/clustagglom/file/names_file.hpp"
#include "cath/clustagglom/get_sorting_scores.hpp"
#include "cath/clustagglom/hierarchy.hpp"
#include "cath/clustagglom/make_clusters_from_merges.hpp"
#include "cath/clustagglom/merge.hpp"
#include "cath/common/container/id_of_str_bidirnl.hpp"
#include "cath/common/file/temp_file.hpp"
#include "cath/test/predicate/files_equal.hpp"

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;

using ::std::filesystem::path;

namespace {

		/// \brief The make_clusters_from_merges_test_suite_fixture to assist in testing make_clusters_from_merges
		struct make_clusters_from_merges_test_suite_fixture : protected clustagglom_fixture {
		protected:
			~make_clusters_from_merges_test_suite_fixture() noexcept = default;

			/// \brief The temp file
			temp_file temp_clusters{ ".cath_tools_test_temp_file.make_clusters_from_merges.%%%%.%%%%-%%%%-%%%%-%%%%" };

			/// \brief The path of the temp_file
			const path temp_clusters_file{ get_filename( temp_clusters ) };

			/// \brief Test that making clusters from the data in the specified names/merges files at
			///        the specified cutoffs generates data matching that in the specified expected file
			void test_clusters_of_merges(const path         &prm_names_file,   ///< The file containing the names data
			                             const path         &prm_merge_file,   ///< The file containing the merges data
			                             const strength_vec &prm_cutoffs,      ///< The cutoffs at which the clusters should be performed
			                             const path         &prm_expected_file ///< The expected file to compare against
			                             ) const {
				const merge_vec          merges                = read_merge_list( prm_merge_file );
				const auto               parse_names_results   = parse_names( prm_names_file );
				const id_of_str_bidirnl &the_id_of_str_bidirnl = parse_names_results.second;
				const size_vec           sorting_indices       = get_sorting_scores(
					the_id_of_str_bidirnl,
					parse_names_results.first
				);

				const hierarchy the_hierarchy = make_clusters_from_merges_and_sort(
					merges,
					sorting_indices,
					prm_cutoffs
				);

				write_cluster( temp_clusters_file, the_hierarchy, the_id_of_str_bidirnl );

				BOOST_CHECK_FILES_EQUAL( temp_clusters_file, prm_expected_file );
			}
		};
} // namespace

BOOST_FIXTURE_TEST_SUITE(make_clusters_from_merges_test_suite, make_clusters_from_merges_test_suite_fixture)

BOOST_AUTO_TEST_CASE(standard_examples_make_clusters_from_merges_correctly) {
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.10.8.260.names",
		CLUSTAGGLOM_DIR() / "1.10.8.260.expected_mergelist",
		{ -100.0, -95.0, -60.0, -35.0 },
		CLUSTAGGLOM_DIR() / "1.10.8.260.expected_clusters"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.230.names",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.230.expected_mergelist",
		{ -100.0, -95.0, -60.0, -35.0 },
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.230.expected_clusters"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.10.287.230.names",
		CLUSTAGGLOM_DIR() / "1.10.287.230.expected_mergelist",
		{ -100.0, -95.0, -60.0, -35.0 },
		CLUSTAGGLOM_DIR() / "1.10.287.230.expected_clusters"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.1770.names",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.1770.expected_mergelist",
		{ -100.0, -95.0, -60.0, -35.0 },
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.1770.expected_clusters"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.10.287.1770.names",
		CLUSTAGGLOM_DIR() / "1.10.287.1770.expected_mergelist",
		{ -100.0, -95.0, -60.0, -35.0 },
		CLUSTAGGLOM_DIR() / "1.10.287.1770.expected_clusters"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.20.5.1350.names",
		CLUSTAGGLOM_DIR() / "1.20.5.1350.expected_mergelist",
		{ -100.0, -95.0, -60.0, -35.0 },
		CLUSTAGGLOM_DIR() / "1.20.5.1350.expected_clusters"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.20.5.4590.names",
		CLUSTAGGLOM_DIR() / "1.20.5.4590.expected_mergelist",
		{ -100.0, -95.0, -60.0, -35.0 },
		CLUSTAGGLOM_DIR() / "1.20.5.4590.expected_clusters"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.20.1270.130.names",
		CLUSTAGGLOM_DIR() / "1.20.1270.130.expected_mergelist",
		{ -100.0, -95.0, -60.0, -35.0 },
		CLUSTAGGLOM_DIR() / "1.20.1270.130.expected_clusters"
	);
}

BOOST_AUTO_TEST_SUITE_END()
