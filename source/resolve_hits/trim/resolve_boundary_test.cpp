/// \file
/// \brief The resolve_boundary test suite

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

#include "resolve_hits/trim/resolve_boundary.hpp"

using namespace cath::rslv;
using namespace cath::seq;

using std::invalid_argument;

BOOST_AUTO_TEST_SUITE(resolve_boundary_test_suite)

BOOST_AUTO_TEST_CASE(throws_on_misordered_ends) {
	BOOST_CHECK_THROW( resolve_boundary( arrow_after_res( 200 ), 100, arrow_after_res( 100 ), 100 ), invalid_argument );
}

BOOST_AUTO_TEST_CASE(throws_on_non_meeting_boundaries) {
	BOOST_CHECK_THROW( resolve_boundary( arrow_after_res( 100 ),  40, arrow_after_res( 200 ),  40 ), invalid_argument );
}

BOOST_AUTO_TEST_CASE(resolves_half_way_for_simple_symmetric) {
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 100, arrow_after_res( 200 ), 100 ) == arrow_after_res( 150 ), "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(removes_in_ratio_to_trimmed_lengths) {
	static_assert(     resolve_boundary( arrow_after_res( 100 ),  50, arrow_after_res( 200 ),  75 ) == arrow_after_res( 140 ), "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(rounds_sensibly_below_and_above_residue_midpoint) {
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 145, arrow_after_res( 200 ), 155 ) == arrow_after_res( 148 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 146, arrow_after_res( 200 ), 154 ) == arrow_after_res( 149 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 147, arrow_after_res( 200 ), 153 ) == arrow_after_res( 149 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 148, arrow_after_res( 200 ), 152 ) == arrow_after_res( 149 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 149, arrow_after_res( 200 ), 151 ) == arrow_after_res( 150 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 150, arrow_after_res( 200 ), 150 ) == arrow_after_res( 150 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 151, arrow_after_res( 200 ), 149 ) == arrow_after_res( 150 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 152, arrow_after_res( 200 ), 148 ) == arrow_after_res( 151 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 153, arrow_after_res( 200 ), 147 ) == arrow_after_res( 151 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 154, arrow_after_res( 200 ), 146 ) == arrow_after_res( 151 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 155, arrow_after_res( 200 ), 145 ) == arrow_after_res( 152 ), "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(gives_spare_residue_to_shorter_trim_at_residue_midpoint) {
	static_assert(     resolve_boundary( arrow_after_res( 100 ),  96, arrow_after_res( 200 ), 104 ) == arrow_after_res( 148 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ),  97, arrow_after_res( 200 ), 103 ) == arrow_after_res( 149 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ),  98, arrow_after_res( 200 ), 102 ) == arrow_after_res( 149 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ),  99, arrow_after_res( 200 ), 101 ) == arrow_after_res( 150 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 100, arrow_after_res( 200 ), 100 ) == arrow_after_res( 150 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 101, arrow_after_res( 200 ),  99 ) == arrow_after_res( 150 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 102, arrow_after_res( 200 ),  98 ) == arrow_after_res( 151 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 103, arrow_after_res( 200 ),  97 ) == arrow_after_res( 151 ), "" );
	static_assert(     resolve_boundary( arrow_after_res( 100 ), 104, arrow_after_res( 200 ),  96 ) == arrow_after_res( 152 ), "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(rounds_up_at_midpoint_with_equal_trims) {
	static_assert( resolve_boundary( arrow_after_res( 100 ), 100, arrow_after_res( 199 ), 100 ) == arrow_after_res( 150 ), "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
