/// \file
/// \brief The regions_limiter test suite

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

#include "biocore/residue_id.hpp"
#include "chopping/region/region.hpp"
#include "chopping/region/regions_limiter.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common;

BOOST_AUTO_TEST_SUITE(regions_limiter_test_suite)


BOOST_AUTO_TEST_CASE(ctor_throws_given_region_without_chain_label) {
	const region_vec the_regions = { region{ 1, 2 } };
	BOOST_CHECK_THROW( regions_limiter the_limiter{ the_regions }, invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(empty) {
	regions_limiter the_limiter{};
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 1 ) ), true );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 2 ) ), true );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 1 ) ), true );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 2 ) ), true );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'C', 1 ) ), true );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'C', 2 ) ), true );
}


BOOST_AUTO_TEST_CASE(whole_chain_segment) {
	const region_vec the_regions = { make_simple_region( 'B' ) };
	regions_limiter the_limiter{ the_regions };
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 1 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 2 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 1 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 2 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'C', 1 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'C', 2 ) ), false );
}

BOOST_AUTO_TEST_CASE(consecutive_chains) {
	const region_vec the_regions = { make_simple_region( 'B' ), make_simple_region( 'C' ) };
	regions_limiter the_limiter{ the_regions };
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 1 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 2 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 1 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 2 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'C', 1 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'C', 2 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'D', 1 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'D', 2 ) ), false );
}

BOOST_AUTO_TEST_CASE(consecutive_segments) {
	const region_vec the_regions = { make_simple_region( 'A', 2, 3 ), make_simple_region( 'A', 4, 5 ) };
	regions_limiter the_limiter{ the_regions };
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 1 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 2 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 3 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 4 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 5 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 6 ) ), false );
}

BOOST_AUTO_TEST_CASE(multiple_segments) {
	const region_vec the_regions = { make_simple_region( 'B', 2, 4 ), make_simple_region( 'C', 1, 1 ) };
	regions_limiter the_limiter{ the_regions };
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 1 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'A', 2 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 1 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 2 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 3 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 4 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'B', 5 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'C', 1 ) ), true  );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'C', 2 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'D', 1 ) ), false );
	BOOST_CHECK_EQUAL( the_limiter.update_residue_is_included( make_residue_id( 'D', 2 ) ), false );
}

BOOST_AUTO_TEST_SUITE_END()
