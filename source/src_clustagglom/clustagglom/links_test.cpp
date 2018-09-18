/// \file
/// \brief The links test suite

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

#include <boost/test/auto_unit_test.hpp>

#include "clustagglom/clustagglom_fixture.hpp"
#include "clustagglom/file/dissimilarities_file.hpp"
#include "clustagglom/file/names_file.hpp"
#include "clustagglom/links.hpp"
#include "common/container/id_of_str_bidirnl.hpp"
#include "common/file/temp_file.hpp"
#include "test/predicate/files_equal.hpp"
#include "clustagglom/get_sorting_scores.hpp"
// #include "clustagglom/hierarchy.hpp"
// #include "clustagglom/make_clusters_from_merges.hpp"
// #include "clustagglom/merge.hpp"

namespace cath { namespace test { } }

// using namespace cath;
using namespace cath::clust;
using namespace cath::common;
using namespace cath::test;

using boost::filesystem::path;

namespace cath {
	namespace test {

		/// \brief The links_test_suite_test_suite_fixture to assist in testing links_test_suite
		struct links_test_suite_test_suite_fixture : protected clustagglom_fixture {
		protected:
			~links_test_suite_test_suite_fixture() noexcept = default;

			/// \brief The temp file
			temp_file temp_links{ ".cath_tools_test_temp_file.links_test_suite.%%%%.%%%%-%%%%-%%%%-%%%%" };

			/// \brief The path of the temp_file
			const path temp_links_file{ get_filename( temp_links ) };

			/// \brief Test that making clusters from the data in the specified names/merges files at
			///        the specified cutoffs generates data matching that in the specified expected file
			void test_clusters_of_merges(const path         &prm_names_file,   ///< The file containing the names data
			                             const path         &prm_links_file,   ///< The file containing the links data
			                             const path         &prm_expected_file ///< The expected file to compare against
			                             ) const {
				id_of_str_bidirnl id_namer;

				const doub_vec props        = parse_names( prm_names_file, id_namer );
				auto           the_links    = parse_dissimilarities( prm_links_file, id_namer, link_dirn::STRENGTH );
				const size_vec sort_indices = get_sorting_scores( id_namer, props );

				write_ordered_links( temp_links_file, the_links, id_namer, sort_indices );

				BOOST_CHECK_FILES_EQUAL( temp_links_file, prm_expected_file );
			}
		};
	}
}

BOOST_FIXTURE_TEST_SUITE(links_test_suite, links_test_suite_test_suite_fixture)

BOOST_AUTO_TEST_CASE(standard_examples_reorder_links_correctly) {
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.10.8.260.names",
		CLUSTAGGLOM_DIR() / "1.10.8.260.nwresults",
		CLUSTAGGLOM_DIR() / "1.10.8.260.expected_reordered_nwresults"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.230.names",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.230.nwresults",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.230.expected_reordered_nwresults"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.10.287.230.names",
		CLUSTAGGLOM_DIR() / "1.10.287.230.nwresults",
		CLUSTAGGLOM_DIR() / "1.10.287.230.expected_reordered_nwresults"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.1770.names",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.1770.nwresults",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.1770.expected_reordered_nwresults"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.10.287.1770.names",
		CLUSTAGGLOM_DIR() / "1.10.287.1770.nwresults",
		CLUSTAGGLOM_DIR() / "1.10.287.1770.expected_reordered_nwresults"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.20.5.1350.names",
		CLUSTAGGLOM_DIR() / "1.20.5.1350.nwresults",
		CLUSTAGGLOM_DIR() / "1.20.5.1350.expected_reordered_nwresults"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.20.5.4590.names",
		CLUSTAGGLOM_DIR() / "1.20.5.4590.nwresults",
		CLUSTAGGLOM_DIR() / "1.20.5.4590.expected_reordered_nwresults"
	);
	test_clusters_of_merges(
		CLUSTAGGLOM_DIR() / "1.20.1270.130.names",
		CLUSTAGGLOM_DIR() / "1.20.1270.130.nwresults",
		CLUSTAGGLOM_DIR() / "1.20.1270.130.expected_reordered_nwresults"
	);
}

BOOST_AUTO_TEST_SUITE_END()
