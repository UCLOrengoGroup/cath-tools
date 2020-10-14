/// \file
/// \brief The clust_id_pot test suite

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

#include "cath/clustagglom/detail/clust_id_pot.hpp"

using namespace cath::clust::detail;

BOOST_AUTO_TEST_SUITE(clust_id_pot_test_suite)


BOOST_AUTO_TEST_CASE(const_methods) {
	const clust_id_pot the_clust_pot{ 3 };

	BOOST_TEST( the_clust_pot.has_index                   ( 0 )      );
	BOOST_TEST( the_clust_pot.get_jumbled_nth_index       ( 0 ) >= 0 );
	BOOST_TEST( the_clust_pot.get_jumbled_nth_index       ( 0 ) <= 2 );
	BOOST_TEST( the_clust_pot.get_min_value_excluding_spec( 0 ) <= 2 );
	BOOST_TEST( the_clust_pot.get_min_value_excluding_spec( 1 ) == 0 );
}


BOOST_AUTO_TEST_CASE(basic0) {
	clust_id_pot the_clust_pot{ 3 };

	BOOST_TEST(   the_clust_pot.has_index( 0 ) );
	BOOST_TEST(   the_clust_pot.has_index( 1 ) );
	BOOST_TEST(   the_clust_pot.has_index( 2 ) );
}

BOOST_AUTO_TEST_CASE(basic1) {
	clust_id_pot the_clust_pot{ 3 };

	the_clust_pot.remove_index( 1 );

	BOOST_TEST(   the_clust_pot.has_index( 0 ) );
	BOOST_TEST( ! the_clust_pot.has_index( 1 ) );
	BOOST_TEST(   the_clust_pot.has_index( 2 ) );
}

BOOST_AUTO_TEST_CASE(basic2) {
	clust_id_pot the_clust_pot{ 3 };

	the_clust_pot.remove_index( 2 );

	BOOST_TEST(   the_clust_pot.has_index( 0 ) );
	BOOST_TEST(   the_clust_pot.has_index( 1 ) );
	BOOST_TEST( ! the_clust_pot.has_index( 2 ) );
}


BOOST_AUTO_TEST_SUITE_END()
