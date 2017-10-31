/// \map
/// \brief The overlap_frac_distn test suite

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

#include "cluster/map/overlap_frac_distn.hpp"
#include "test/boost_addenda/boost_check_equal_ranges.hpp"

using namespace cath::clust;
using namespace cath::common;

BOOST_AUTO_TEST_SUITE(overlap_frac_distn_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	overlap_frac_distn bob;

	BOOST_CHECK(   bob.empty() );
	bob.add_overlap_fraction( 0.1 );
	BOOST_CHECK( ! bob.empty() );
	bob.add_overlap_fraction( 0.1 );
	bob.add_overlap_fraction( 0.2 );
	bob.add_overlap_fraction( 0.2 );
	bob.add_overlap_fraction( 0.2 );
	bob.add_overlap_fraction( 0.3 );
	BOOST_CHECK( ! bob.empty() );

	BOOST_CHECK_EQUAL( bob.size(), 6 );
	BOOST_CHECK_EQUAL( bob.get_num_in_range ( 0.00, 0.05 ), 0 );
	BOOST_CHECK_EQUAL( bob.get_num_in_range ( 0.12, 0.15 ), 0 );
	BOOST_CHECK_EQUAL( bob.get_num_in_range ( 0.08, 0.15 ), 2 );
	BOOST_CHECK_EQUAL( bob.get_num_in_range ( 0.08, 0.20 ), 2 );
	BOOST_CHECK_EQUAL( bob.get_num_in_range ( 0.08, 0.21 ), 5 );
	BOOST_CHECK_EQUAL( bob.get_num_in_range ( 0.00, 2.00 ), 6 );

	BOOST_CHECK_EQUAL( bob.get_frac_at_percentile (   0.00 ), 0.1 );
	BOOST_CHECK_EQUAL( bob.get_frac_at_percentile (  25.00 ), 0.1 );
	BOOST_CHECK_EQUAL( bob.get_frac_at_percentile (  50.00 ), 0.2 );
	BOOST_CHECK_EQUAL( bob.get_frac_at_percentile (  75.00 ), 0.2 );
	BOOST_CHECK_EQUAL( bob.get_frac_at_percentile ( 100.00 ), 0.3 );
}

BOOST_AUTO_TEST_CASE(adding) {
	const overlap_frac_distn alice = build_overlap_frac_distn_from_overlap_fractions( { 0.1, 0.1, 0.2, 0.2, 0.2, 0.3 } );
	const overlap_frac_distn bob   = build_overlap_frac_distn_from_overlap_fractions( { 0.6, 0.7, 0.7, 0.8           } );
	const auto charlie = alice + bob;

	BOOST_CHECK_EQUAL( charlie.size(),                            10  );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (   0.00 ), 0.1 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (  10.00 ), 0.1 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (  20.00 ), 0.1 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (  30.00 ), 0.2 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (  40.00 ), 0.2 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (  50.00 ), 0.2 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (  60.00 ), 0.3 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (  70.00 ), 0.6 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (  80.00 ), 0.7 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile (  90.00 ), 0.7 );
	BOOST_CHECK_EQUAL( charlie.get_frac_at_percentile ( 100.00 ), 0.8 );
}

BOOST_AUTO_TEST_SUITE_END()
