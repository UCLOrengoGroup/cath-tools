/// \file
/// \brief The region test suite

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

#include <boost/optional/optional_io.hpp>

#include "chopping/chopping_type_aliases.hpp"
#include "chopping/region/region.hpp"
#include "test/test_tools.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common::test;

using boost::none;

BOOST_AUTO_TEST_SUITE(region_test_suite)

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL( to_string( make_simple_region( 'A', 121, 232 ) ), "region[ chain:A, start_name:121, stop_name:232 ]" );
	BOOST_CHECK_EQUAL( to_string( make_simple_region( 121, 232      ) ), "region[ start_idx:121, stop_idx:232 ]"            );

	BOOST_CHECK_EQUAL( to_string( make_simple_region( 'A'           ) ), "region[ chain:A ]"                                );
}

BOOST_AUTO_TEST_CASE(equality_works) {
	check_equality_operators_on_diff_vals_range( region_vec{
		region{ chain_label{ 'K' } },
		make_simple_region( 121, 232 ),
		make_simple_region( 'K', 121, 232 )
	} );
}

BOOST_AUTO_TEST_CASE(locating_of_whole_chain_region_returns_none) {
	BOOST_TEST( get_residue_locating( make_simple_region( 'A' ) ) == none );
}

BOOST_AUTO_TEST_SUITE_END()
