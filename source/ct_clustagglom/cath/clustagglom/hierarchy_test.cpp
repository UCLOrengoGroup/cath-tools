/// \file
/// \brief The hierarchy test suite

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

#include <boost/test/unit_test.hpp>

#include "cath/clustagglom/hierarchy.hpp"

using namespace ::cath::clust;
using namespace ::cath::common;

using ::boost::test_tools::per_element;

BOOST_AUTO_TEST_SUITE(hierarchy_test_suite)

BOOST_AUTO_TEST_CASE(get_rep_indices_of_single_direct_entry_returns_that_entry) {
	constexpr auto array_of_zero  = make_array( 0_z );
	const hierarchy the_hierarchy = make_hierarchy_from_reversed_without_root_layer(
		{ hierarchy_layer{} },
		array_of_zero
	);
	BOOST_TEST( get_rep_indices( the_hierarchy ) == array_of_zero, per_element{} );
}

BOOST_AUTO_TEST_SUITE_END()
