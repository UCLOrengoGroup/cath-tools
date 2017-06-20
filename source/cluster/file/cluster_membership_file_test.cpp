/// \file
/// \brief The cluster_membership_file test suite

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

#include "cluster/file/cluster_membership_file.hpp"

namespace cath { namespace test { } }

using namespace cath::clust;
using namespace cath::common;
using namespace cath::test;

using std::string;

namespace cath {
	namespace test {

		/// \brief The cluster_membership_file_test_suite_fixture to assist in testing cluster_membership_file
		struct cluster_membership_file_test_suite_fixture {
		protected:
			~cluster_membership_file_test_suite_fixture() noexcept = default;

			/// \brief An example "old" membership string
			const string old_membership_str = R"(369 126c3de59272d9c6d2de490d47f9f19d
523 4f5c5f455b1589578cd42403f846f48b/108-168
523 6b02ac7d3d285bc5100e142c133808ff/66-129
523 14e361063ccad7801d24cff355749dfb/115-178
523 16b30c0d2fe3367d022dbe8dd2f38d2f/87-141
523 16f85e1adbb904c4b56067a7f35b92c4/89-136
523 31f815681481495e0b0f23e5c5937d49/66-129
)";

			/// \brief An example "new" membership string
			const string new_membership_str = R"(335 14e361063ccad7801d24cff355749dfb/206-244
335 72c898c9ba2ec6809a457be7ce0d72fb/162-200
414 4f5c5f455b1589578cd42403f846f48b/200-235
414 3724f6c949f50a3e71dd4166debe1c8e/163-205
416 16f85e1adbb904c4b56067a7f35b92c4/87-135_203-228
416 126c3de59272d9c6d2de490d47f9f19d
)";

			/// \brief A id_of_string for mapping the seq names to ID numbers
			id_of_string seq_id_mapper;

			/// \brief Example old_cluster_data parsed from the old_membership_str
			const old_cluster_data old_data = parse_old_membership( old_membership_str, seq_id_mapper );

			/// \brief Example new_cluster_data parsed from the new_membership_str
			const new_cluster_data new_data = parse_new_membership( new_membership_str, seq_id_mapper );

		};
	}
}

BOOST_FIXTURE_TEST_SUITE(cluster_membership_file_test_suite, cluster_membership_file_test_suite_fixture)

BOOST_AUTO_TEST_CASE(old_cluster_data__to_string) {
	const string expected_str = "old_cluster_data[2 clusters, "
		R"(clusters{ 0("369"): 126c3de59272d9c6d2de490d47f9f19d(*), )"
		R"(1("523"): 4f5c5f455b1589578cd42403f846f48b([108-168]), )"
		"6b02ac7d3d285bc5100e142c133808ff([66-129]), "
		"14e361063ccad7801d24cff355749dfb([115-178]), "
		"16b30c0d2fe3367d022dbe8dd2f38d2f([87-141]), "
		"16f85e1adbb904c4b56067a7f35b92c4([89-136]), "
		"31f815681481495e0b0f23e5c5937d49([66-129]) } ]";

	BOOST_CHECK_EQUAL( to_string( old_data ), expected_str );
}

BOOST_AUTO_TEST_CASE(new_cluster_data__to_string) {
	const string expected_str = "new_cluster_data[3 clusters, "
		R"(cluster_sizes{ 0("335"):2, 1("414"):2, 2("416"):2 }, )"
		  "cluster_index_by_seq_regions{ "
		  "126c3de59272d9c6d2de490d47f9f19d:(*->2), "
		  "14e361063ccad7801d24cff355749dfb:([206-244]->0), "
		  "16f85e1adbb904c4b56067a7f35b92c4:([87-135,203-228]->2), "
		  "3724f6c949f50a3e71dd4166debe1c8e:([163-205]->1), "
		  "4f5c5f455b1589578cd42403f846f48b:([200-235]->1), "
		  "72c898c9ba2ec6809a457be7ce0d72fb:([162-200]->0) } ]";

	BOOST_CHECK_EQUAL( to_string( new_data ), expected_str );
}

BOOST_AUTO_TEST_SUITE_END()
