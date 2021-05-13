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

#include "cath/biocore/residue_id.hpp"
#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/test/boost_test_print_type.hpp"
#include "cath/test/test_tools.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common::test;

using ::std::nullopt;

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
	BOOST_TEST( get_residue_locating( make_simple_region( 'A' ) ) == nullopt );
}

BOOST_AUTO_TEST_CASE( get_start_id_and_get_stop_id_work ) {
	BOOST_TEST( get_start_id( make_simple_region( 'K', 121, 232 ) ) == make_residue_id( 'K', 121 ) );
	BOOST_TEST( get_stop_id( make_simple_region( 'K', 121, 232 ) ) == make_residue_id( 'K', 232 ) );
}

BOOST_AUTO_TEST_SUITE_END()
