/// \file
/// \brief The calc_complete_linkage_merge_list test suite

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

#include "clustagglom/calc_complete_linkage_merge_list.hpp"
#include "clustagglom/clustagglom_fixture.hpp"
#include "clustagglom/file/dissimilarities_file.hpp"
#include "clustagglom/file/names_file.hpp"
#include "clustagglom/get_sorting_scores.hpp"
#include "clustagglom/links.hpp"
#include "clustagglom/merge.hpp"
#include "common/container/id_of_str_bidirnl.hpp"
#include "common/file/temp_file.hpp"
#include "test/predicate/files_equal.hpp"

using namespace cath::clust;
using namespace cath::common;
using namespace cath::test;

using boost::filesystem::path;
using boost::test_tools::per_element;
using std::move;
using std::numeric_limits;

namespace cath {
	namespace test {

		/// \brief The calc_complete_linkage_merge_list_test_suite_fixture to assist in testing calc_complete_linkage_merge_list
		struct calc_complete_linkage_merge_list_test_suite_fixture : protected clustagglom_fixture {
		protected:
			~calc_complete_linkage_merge_list_test_suite_fixture() noexcept = default;

			/// \brief The temp file
			temp_file temp_mergelist{ ".cath_tools_test_temp_file.calc_complete_linkage_merge_list.%%%%.%%%%-%%%%-%%%%-%%%%" };

			/// \brief The path of the temp_file
			const path temp_mergelist_file{ get_filename( temp_mergelist ) };

			/// \brief Test that making the complete-linkage merge list for the data from the specified names/links files
			///        (up to the specified max_dissim and reading the links from the file in the specified direction)
			///        generates  data matching that in the specified expected file
			void test_complete_linkage_merge_list(const path      &prm_names_file,                                             ///< The file containing the names data
			                                      const path      &prm_links_file,                                             ///< The file containing the links data
			                                      const path      &prm_expected_file,                                          ///< The expected file to compare against
			                                      const strength  &prm_max_dissim = std::numeric_limits<strength>::infinity(), ///< The maximum dissimilarity at which merges may still happen
			                                      const link_dirn &prm_link_dirn  = link_dirn::STRENGTH                        ///< Whether the links in the input file represent strengths or dissimilarities
			                                      ) const {
				id_of_str_bidirnl the_id_of_str_bidirnl;
				const auto        props           = parse_names( prm_names_file, the_id_of_str_bidirnl );
				auto              dissims         = parse_dissimilarities( prm_links_file, the_id_of_str_bidirnl, prm_link_dirn );
				const size_vec    sorting_indices = get_sorting_scores( the_id_of_str_bidirnl, props );
				const auto        the_merge_list  = calc_complete_linkage_merge_list(
					move( dissims ),
					sorting_indices,
					prm_max_dissim
				);

				write_merge_list( temp_mergelist_file, the_merge_list );

				BOOST_CHECK_FILES_EQUAL( temp_mergelist_file, prm_expected_file );
			}

		};
	}
}

BOOST_FIXTURE_TEST_SUITE(calc_complete_linkage_merge_list_test_suite, calc_complete_linkage_merge_list_test_suite_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	const links input_links = make_links( item_item_strength_tpl_vec{ {
		item_item_strength_tpl{ 0, 1, 1.0 },
		item_item_strength_tpl{ 1, 2, 1.0 },
		item_item_strength_tpl{ 2, 3, 1.0 },
	} } );

	const merge_vec expected = merge_vec{ {
		merge{ 0, 1, 4, 1.0 },
		merge{ 2, 3, 5, 1.0 },
	} };

	BOOST_TEST( calc_complete_linkage_merge_list( input_links, 4, 3.0 ) == expected, per_element{} );
}

BOOST_AUTO_TEST_CASE(standard_examples_complete_linkage_merge_correctly) {
	test_complete_linkage_merge_list(
		CLUSTAGGLOM_DIR() / "1.10.8.260.names",
		CLUSTAGGLOM_DIR() / "1.10.8.260.nwresults",
		CLUSTAGGLOM_DIR() / "1.10.8.260.expected_mergelist",
		35.0,
		link_dirn::STRENGTH
	);
	test_complete_linkage_merge_list(
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.230.names",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.230.nwresults",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.230.expected_mergelist",
		35.0,
		link_dirn::STRENGTH
	);
	test_complete_linkage_merge_list(
		CLUSTAGGLOM_DIR() / "1.10.287.230.names",
		CLUSTAGGLOM_DIR() / "1.10.287.230.nwresults",
		CLUSTAGGLOM_DIR() / "1.10.287.230.expected_mergelist",
		35.0,
		link_dirn::STRENGTH
	);
	test_complete_linkage_merge_list(
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.1770.names",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.1770.nwresults",
		CLUSTAGGLOM_DIR() / "simplified_1.10.287.1770.expected_mergelist",
		35.0,
		link_dirn::STRENGTH
	);
	test_complete_linkage_merge_list(
		CLUSTAGGLOM_DIR() / "1.10.287.1770.names",
		CLUSTAGGLOM_DIR() / "1.10.287.1770.nwresults",
		CLUSTAGGLOM_DIR() / "1.10.287.1770.expected_mergelist",
		35.0,
		link_dirn::STRENGTH
	);
	test_complete_linkage_merge_list(
		CLUSTAGGLOM_DIR() / "1.20.5.1350.names",
		CLUSTAGGLOM_DIR() / "1.20.5.1350.nwresults",
		CLUSTAGGLOM_DIR() / "1.20.5.1350.expected_mergelist",
		35.0,
		link_dirn::STRENGTH
	);
	test_complete_linkage_merge_list(
		CLUSTAGGLOM_DIR() / "1.20.5.4590.names",
		CLUSTAGGLOM_DIR() / "1.20.5.4590.nwresults",
		CLUSTAGGLOM_DIR() / "1.20.5.4590.expected_mergelist",
		35.0,
		link_dirn::STRENGTH
	);
	test_complete_linkage_merge_list(
		CLUSTAGGLOM_DIR() / "1.20.1270.130.names",
		CLUSTAGGLOM_DIR() / "1.20.1270.130.nwresults",
		CLUSTAGGLOM_DIR() / "1.20.1270.130.expected_mergelist",
		35.0,
		link_dirn::STRENGTH
	);
}

BOOST_AUTO_TEST_SUITE_END()
