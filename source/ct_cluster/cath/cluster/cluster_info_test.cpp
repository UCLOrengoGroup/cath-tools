/// \file
/// \brief The cluster_info test suite

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

#include <boost/algorithm/cxx11/is_sorted.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/cluster/cluster_info.hpp"

namespace cath { namespace test { } }

using namespace cath::clust;
using namespace cath::common;
using namespace cath::seq;
using namespace cath::test;

using boost::algorithm::is_strictly_increasing;
using boost::none;
using std::pair;
using std::string;
using std::vector;

namespace cath {
	namespace test {

		/// \brief The cluster_info_test_suite_fixture to assist in testing cluster_info
		struct cluster_info_test_suite_fixture {
		protected:
			~cluster_info_test_suite_fixture() noexcept = default;

			/// \brief Type alias for a pair of string and seq_seg_run_opt
			using str_seq_seg_run_opt_pair     = pair  < string, seq_seg_run_opt  >;

			/// \brief Type alias for a vector of str_seq_seg_run_opt_pair
			using str_seq_seg_run_opt_pair_vec = vector< str_seq_seg_run_opt_pair >;

			/// \brief Make a cluster_info from the specified str_seq_seg_run_opt_pair_vec
			cluster_info make_cluster_info(const str_seq_seg_run_opt_pair_vec &prm_data ///< The data from which to build the cluster_info
			                               ) {
				cluster_info result{};

				// \TODO Come C++17 and structured bindings, use here
				for (const pair<string, seq_seg_run_opt> &datum : prm_data) {
					result.add_entry( datum.first, datum.second );
				}
				return result;
			}
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(cluster_info_test_suite, cluster_info_test_suite_fixture)


BOOST_AUTO_TEST_CASE(has_correct_properties_after_construction) {
	const cluster_info cluster_info_a{};
	BOOST_CHECK_EQUAL( cluster_info_a.get_size                 (),    0                          );
	BOOST_CHECK_THROW( cluster_info_a.get_total_sqrt_length    (),    invalid_argument_exception );
	BOOST_CHECK_THROW( cluster_info_a.get_total_mid_point_index(),    invalid_argument_exception );
	BOOST_CHECK_THROW( cluster_info_a.get_lowest_domain_id     (),    invalid_argument_exception );
	BOOST_CHECK_THROW( get_average_mid_point_index( cluster_info_a ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(has_correct_properties_after_adding_two_with_no_segs) {
	const cluster_info cluster_info_a = make_cluster_info( {
		str_seq_seg_run_opt_pair{ "sue",  none },
		str_seq_seg_run_opt_pair{ "mary", none },
	} );
	BOOST_CHECK_EQUAL( cluster_info_a.get_size                 (),       2 );
	BOOST_CHECK_EQUAL( cluster_info_a.get_total_sqrt_length    (),    none );
	BOOST_CHECK_EQUAL( cluster_info_a.get_total_mid_point_index(),    none );
	BOOST_CHECK_EQUAL( cluster_info_a.get_lowest_domain_id     (),  "mary" );
	BOOST_CHECK_EQUAL( get_average_mid_point_index( cluster_info_a ), none );
}

BOOST_AUTO_TEST_CASE(has_correct_properties_after_adding_two_with_segs) {
	const cluster_info cluster_info_a = make_cluster_info( {
		str_seq_seg_run_opt_pair{ "sue",  make_seq_seg_run_from_res_indices(  10,  90 ) },
		str_seq_seg_run_opt_pair{ "mary", make_seq_seg_run_from_res_indices( 110, 190 ) },
	} );
	BOOST_CHECK_EQUAL( cluster_info_a.get_size                    (),      2 );
	BOOST_CHECK_EQUAL( cluster_info_a.get_total_sqrt_length       (),   18.0 );
	BOOST_CHECK_EQUAL( cluster_info_a.get_total_mid_point_index   (),  200.0 );
	BOOST_CHECK_EQUAL( cluster_info_a.get_lowest_domain_id        (), "mary" );
	BOOST_CHECK_EQUAL( get_average_mid_point_index( cluster_info_a ),  100.0 );
}

BOOST_AUTO_TEST_CASE(treats_empty_string_as_lower) {
	const cluster_info cluster_info_sue_then_empty = make_cluster_info( {
		str_seq_seg_run_opt_pair{ "sue",  none },
		str_seq_seg_run_opt_pair{ "",     none },
	} );
	BOOST_CHECK_EQUAL( cluster_info_sue_then_empty.get_lowest_domain_id(), "" );
	const cluster_info cluster_info_empty_then_sue = make_cluster_info( {
		str_seq_seg_run_opt_pair{ "",     none },
		str_seq_seg_run_opt_pair{ "sue",  none },
	} );
	BOOST_CHECK_EQUAL( cluster_info_empty_then_sue.get_lowest_domain_id(), "" );
}

BOOST_AUTO_TEST_CASE(less_than_is_correct) {
	const vector<cluster_info> cluster_infos = {

		make_cluster_info( {
			str_seq_seg_run_opt_pair{ "sue", make_seq_seg_run_from_res_indices(  10, 990 ) },
		} ),


		make_cluster_info( {
			str_seq_seg_run_opt_pair{ "rae", make_seq_seg_run_from_res_indices(  10,  90 ) },
			str_seq_seg_run_opt_pair{ "que", make_seq_seg_run_from_res_indices(  10,  90 ) },
		} ),

		make_cluster_info( {
			str_seq_seg_run_opt_pair{ "pam", make_seq_seg_run_from_res_indices(  10, 333 ) },
		} ),

		make_cluster_info( {
			str_seq_seg_run_opt_pair{ "ona", make_seq_seg_run_from_res_indices( 110, 433 ) },
		} ),

		make_cluster_info( {
			str_seq_seg_run_opt_pair{ "onb", make_seq_seg_run_from_res_indices( 110, 433 ) },
		} ),

		make_cluster_info( {
			str_seq_seg_run_opt_pair{ "nas",  none },
			str_seq_seg_run_opt_pair{ "may",  none },
		} ),

		make_cluster_info( {
			str_seq_seg_run_opt_pair{ "kim",  none },
		} ),

		make_cluster_info( {
			str_seq_seg_run_opt_pair{ "liz",  none },
		} ),

	};
	BOOST_CHECK( is_strictly_increasing( cluster_infos ) );
}


BOOST_AUTO_TEST_SUITE_END()
