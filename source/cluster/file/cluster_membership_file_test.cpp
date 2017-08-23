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

#include <boost/optional/optional_io.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "cluster/file/cluster_membership_file.hpp"
#include "cluster/map/map_clusters.hpp"
#include "cluster/options/spec/clust_mapping_spec.hpp"
#include "common/boost_addenda/test/boost_check_no_throw_diag.hpp"
#include "exception/runtime_error_exception.hpp"

namespace cath { namespace test { } }

using namespace cath::clust;
using namespace cath::common;
using namespace cath::test;

using boost::none;
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
414 4f5c5f455b1589578cd42403f846f48b/110-173
414 3724f6c949f50a3e71dd4166debe1c8e/163-205
416 3457216f85e1adbb904c4b56067a2457
416 126c3de59272d9c6d2de490d47f9f19d
)";

			/// \brief A id_of_string for mapping the seq names to ID numbers
			// cluster_name_ider seq_id_mapper; remove include of id_of_string
			id_of_string seq_id_mapper;

			/// \brief Example old_cluster_data parsed from the old_membership_str
			const old_cluster_data old_data = parse_old_membership( old_membership_str, seq_id_mapper );

			/// \brief Example new_cluster_data parsed from the new_membership_str
			const new_cluster_data new_data = parse_new_membership( new_membership_str, seq_id_mapper );

			/// \brief Example of input that's invalid due to only having one entry per line
			const string single_column_input_str         = "a\nb\nc\n";

			/// \brief Example of input that's invalid due to only having one entry per line (even though it's surrounded by space characters)
			const string space_single_column_input_str   = " a \n b \n c \n";

			/// \brief Example of input that's valid despite having trailing whitespace
			const string trailing_space_input_str        = "a 1 \nb 2 \nc 3 \n";

			/// \brief Example of input that's invalid due to having a spurious extra column
			const string spurious_extra_column_input_str = "a 1 x\nb 2 y\nc 3 z\n";

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

BOOST_AUTO_TEST_CASE(old_cluster_info_first) {
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( old_data, 0 ).get_size                    (),                                          1 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( old_data, 0 ).get_total_sqrt_length       (),                                       none );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( old_data, 0 ).get_total_mid_point_index   (),                                       none );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( old_data, 0 ).get_lowest_domain_id        (),         "126c3de59272d9c6d2de490d47f9f19d" );
	BOOST_CHECK_EQUAL( get_average_mid_point_index( get_info_of_cluster_of_id( old_data, 0 ) ),                                       none );
}

BOOST_AUTO_TEST_CASE(old_cluster_info_second) {
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( old_data, 1 ).get_size                    (),                                          6 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( old_data, 1 ).get_total_sqrt_length       (),                         46.154651393277824 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( old_data, 1 ).get_total_mid_point_index   (),                                      706.0 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( old_data, 1 ).get_lowest_domain_id        (), "14e361063ccad7801d24cff355749dfb/115-178" );
	BOOST_CHECK_EQUAL( get_average_mid_point_index( get_info_of_cluster_of_id( old_data, 1 ) ),                         117.66666666666667 );
}

BOOST_AUTO_TEST_CASE(new_cluster_data__to_string) {
	const string expected_str = "new_cluster_data[3 clusters, "
		R"(cluster_sizes{ 0("335"):2, 1("414"):2, 2("416"):2 }, )"
		  "cluster_index_by_seq_regions{ "
		  "126c3de59272d9c6d2de490d47f9f19d:(*->2), "
		  "14e361063ccad7801d24cff355749dfb:([206-244]->0), "
		  "3457216f85e1adbb904c4b56067a2457:(*->2), "
		  "3724f6c949f50a3e71dd4166debe1c8e:([163-205]->1), "
		  "4f5c5f455b1589578cd42403f846f48b:([110-173]->1), "
		  "72c898c9ba2ec6809a457be7ce0d72fb:([162-200]->0) } ]";

	BOOST_CHECK_EQUAL( to_string( new_data ), expected_str );
}

BOOST_AUTO_TEST_CASE(new_cluster_info_first) {
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 0 ).get_size                    (),                                          2 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 0 ).get_total_sqrt_length       (),                         12.489995996796796 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 0 ).get_total_mid_point_index   (),                                      406.0 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 0 ).get_lowest_domain_id        (), "14e361063ccad7801d24cff355749dfb/206-244" );
	BOOST_CHECK_EQUAL( get_average_mid_point_index( get_info_of_cluster_of_id( new_data, 0 ) ),                                      203.0 );
}

BOOST_AUTO_TEST_CASE(new_cluster_info_second) {
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 1 ).get_size                    (),                                          2 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 1 ).get_total_sqrt_length       (),                         14.557438524302000 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 1 ).get_total_mid_point_index   (),                                      325.5 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 1 ).get_lowest_domain_id        (), "3724f6c949f50a3e71dd4166debe1c8e/163-205" );
	BOOST_CHECK_EQUAL( get_average_mid_point_index( get_info_of_cluster_of_id( new_data, 1 ) ),                                     162.75 );
}

BOOST_AUTO_TEST_CASE(new_cluster_info_third) {
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 2 ).get_size                    (),                                          2 );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 2 ).get_total_sqrt_length       (),                                       none );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 2 ).get_total_mid_point_index   (),                                       none );
	BOOST_CHECK_EQUAL( get_info_of_cluster_of_id( new_data, 2 ).get_lowest_domain_id        (),         "126c3de59272d9c6d2de490d47f9f19d" );
	BOOST_CHECK_EQUAL( get_average_mid_point_index( get_info_of_cluster_of_id( new_data, 2 ) ),                                       none );
}

// BOOST_AUTO_TEST_CASE(map_clusters_does_stuff) {
// 	map_clusters( old_data, new_data, clust_mapping_spec{} );
// 	BOOST_CHECK( true );
// }


BOOST_AUTO_TEST_SUITE(edge_case_input)



BOOST_AUTO_TEST_SUITE(throws_on_single_column)

BOOST_AUTO_TEST_CASE(throws_on_old_parse_of_single_column) {
	BOOST_CHECK_THROW( parse_old_membership( single_column_input_str, seq_id_mapper ), runtime_error_exception );
}

BOOST_AUTO_TEST_CASE(throws_on_new_parse_of_single_column) {
	BOOST_CHECK_THROW( parse_new_membership( single_column_input_str, seq_id_mapper ), runtime_error_exception );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(throws_on_space_single_column)

BOOST_AUTO_TEST_CASE(throws_on_old_parse_of_space_single_column) {
	BOOST_CHECK_THROW( parse_old_membership( space_single_column_input_str, seq_id_mapper ), runtime_error_exception );
}

BOOST_AUTO_TEST_CASE(throws_on_new_parse_of_space_single_column) {
	BOOST_CHECK_THROW( parse_new_membership( space_single_column_input_str, seq_id_mapper ), runtime_error_exception );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(accepts_trailing_space)

BOOST_AUTO_TEST_CASE(old_parse_accepts_trailing_space) {
	BOOST_CHECK_NO_THROW_DIAG( parse_old_membership( trailing_space_input_str, seq_id_mapper ) );
}

BOOST_AUTO_TEST_CASE(new_parse_accepts_trailing_space) {
	BOOST_CHECK_NO_THROW_DIAG( parse_new_membership( trailing_space_input_str, seq_id_mapper ) );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(rejects_spurious_extra_column)

BOOST_AUTO_TEST_CASE(throws_on_old_parse_of_spurious_extra_column) {
	BOOST_CHECK_THROW( parse_old_membership( spurious_extra_column_input_str, seq_id_mapper ), runtime_error_exception );
}

BOOST_AUTO_TEST_CASE(throws_on_new_parse_of_spurious_extra_column) {
	BOOST_CHECK_THROW( parse_new_membership( spurious_extra_column_input_str, seq_id_mapper ), runtime_error_exception );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()
