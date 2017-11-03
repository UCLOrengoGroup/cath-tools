/// \map
/// \brief The map_results test suite

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

#include "cluster/map/map_results.hpp"
#include "common/size_t_literal.hpp"
#include "common/type_aliases.hpp"
#include "common/exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::clust;
using namespace cath::common;

BOOST_AUTO_TEST_SUITE(map_results_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( get_num_mapped_entries( map_results{ {}, {}, {}, {          }, {          }, 0, {}, {}, {} } ), 0                          );
	BOOST_CHECK_EQUAL( get_num_mapped_entries( map_results{ {}, {}, {}, { 5_z, 4_z }, { 3_z, 6_z }, 0, {}, {}, {} } ), 9                          );
	BOOST_CHECK_THROW( get_num_mapped_entries( map_results{ {}, {}, {}, { 5_z, 4_z }, { 3_z, 7_z }, 0, {}, {}, {} } ), invalid_argument_exception );
}

// BOOST_AUTO_TEST_CASE(basic) {
// 	const map_results the_results {
// 		potential_map_vec  {} /* chosen_maps                  */ ,
// 		potential_map_vec  {} /* other_maps                   */ ,
// 		size_vec           {} /* unmapped_new_cluster_indices */ ,
// 		size_vec           {}  num_mapped_by_new_cluster     ,
// 		size_vec           {} /* num_mapped_by_old_cluster    */ ,
// 		doub_vec           {} /* domain_mapping_fractions     */ ,
// 		doub_vec           {} /* cluster_mapping_fractions    */ ,
// 		clust_mapping_spec {} /* the_spec                     */
// 	};
// }

BOOST_AUTO_TEST_SUITE_END()
