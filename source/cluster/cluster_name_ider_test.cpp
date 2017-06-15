/// \file
/// \brief The cluster_name_ider test suite

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

#include "cluster/cluster_name_ider.hpp"

using namespace cath::clust;

using boost::make_optional;
using boost::none;

BOOST_AUTO_TEST_SUITE(cluster_name_ider_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	cluster_name_ider the_ider;

	BOOST_CHECK      (   the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 0 );

	BOOST_CHECK_EQUAL(   the_ider.add_name      ( "motorcycle" ), 0            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 0            ), "motorcycle" );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "motorcycle" ), 0            );
	BOOST_CHECK      ( ! the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 1 );

	BOOST_CHECK_EQUAL(   the_ider.add_name      ( "emptiness"  ), 1            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 0            ), "motorcycle" );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "motorcycle" ), 0            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 1            ), "emptiness"  );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "emptiness"  ), 1            );
	BOOST_CHECK      ( ! the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 2 );

	BOOST_CHECK_EQUAL(   the_ider.add_name      ( "motorcycle" ), 0            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 0            ), "motorcycle" );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "motorcycle" ), 0            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 1            ), "emptiness"  );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "emptiness"  ), 1            );
	BOOST_CHECK      ( ! the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 2 );

	the_ider.clear();

	BOOST_CHECK      (   the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 0 );
}


BOOST_AUTO_TEST_SUITE(largest_number_if_names_all_numeric_integers_fn)

BOOST_AUTO_TEST_CASE(returns_largest_negative) {
	cluster_name_ider the_ider;

	the_ider.add_name( "-6"  );
	the_ider.add_name( "-7"  );
	the_ider.add_name( "-9"  );
	the_ider.add_name( "-8"  );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), make_optional( static_cast<ptrdiff_t>( -6 ) ) );
}

BOOST_AUTO_TEST_CASE(handles_non_numeric) {
	cluster_name_ider the_ider;

	the_ider.add_name( "-6"  );
	the_ider.add_name( "-7"  );
	the_ider.add_name( "-9"  );
	the_ider.add_name( "-8"  );
	the_ider.add_name( "bob" );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), none );
}


BOOST_AUTO_TEST_CASE(rejects_scientific_notation) {
	cluster_name_ider the_ider;

	the_ider.add_name( "1e2"  );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), none );
}

BOOST_AUTO_TEST_CASE(rejects_empty) {
	cluster_name_ider the_ider;

	the_ider.add_name( ""  );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), none );
}

BOOST_AUTO_TEST_CASE(rejects_single_dash) {
	cluster_name_ider the_ider;

	the_ider.add_name( "-"  );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), none );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()
