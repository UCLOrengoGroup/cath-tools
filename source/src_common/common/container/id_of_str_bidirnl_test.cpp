/// \file
/// \brief The id_of_str_bidirnl test suite

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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
#include <boost/test/unit_test.hpp>

#include "common/container/id_of_str_bidirnl.hpp"

#include <string>

using namespace cath::common;
using namespace std::literals::string_literals;

using boost::make_optional;
using boost::none;
using boost::string_ref;
using std::string;

BOOST_AUTO_TEST_SUITE(id_of_str_bidirnl_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	id_of_str_bidirnl the_ider;

	BOOST_CHECK      (   the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 0 );

	BOOST_CHECK_EQUAL(   the_ider.add_name      ( "motorcycle"s ), 0             );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 0             ), "motorcycle"s );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "motorcycle"s ), 0             );
	BOOST_CHECK      ( ! the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 1 );

	BOOST_CHECK_EQUAL(   the_ider.add_name      ( "emptiness"s  ), 1             );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 0             ), "motorcycle"s );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "motorcycle"s ), 0             );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 1             ), "emptiness"s  );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "emptiness"s  ), 1             );
	BOOST_CHECK      ( ! the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 2 );

	BOOST_CHECK_EQUAL(   the_ider.add_name      ( "motorcycle"s ), 0             );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 0             ), "motorcycle"s );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "motorcycle"s ), 0             );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 1             ), "emptiness"s  );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( "emptiness"s  ), 1             );
	BOOST_CHECK      ( ! the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 2 );

	the_ider.clear();

	BOOST_CHECK      (   the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 0 );
}

BOOST_AUTO_TEST_CASE(from_string_views) {
	const string motorcycle_str{ "motorcycle" };
	const string emptiness_str { "emptiness"  };

	id_of_str_bidirnl the_ider;

	BOOST_CHECK      (   the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 0 );

	BOOST_CHECK_EQUAL(   the_ider.add_name      ( string_ref{ motorcycle_str } ), 0                            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 0                            ), string_ref{ motorcycle_str } );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( string_ref{ motorcycle_str } ), 0                            );
	BOOST_CHECK      ( ! the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 1 );

	BOOST_CHECK_EQUAL(   the_ider.add_name      ( string_ref{ emptiness_str  } ), 1                            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 0                            ), string_ref{ motorcycle_str } );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( string_ref{ motorcycle_str } ), 0                            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 1                            ), string_ref{ emptiness_str  } );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( string_ref{ emptiness_str  } ), 1                            );
	BOOST_CHECK      ( ! the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 2 );

	BOOST_CHECK_EQUAL(   the_ider.add_name      ( string_ref{ motorcycle_str } ), 0                            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 0                            ), string_ref{ motorcycle_str } );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( string_ref{ motorcycle_str } ), 0                            );
	BOOST_CHECK_EQUAL(   the_ider.get_name_of_id( 1                            ), string_ref{ emptiness_str  } );
	BOOST_CHECK_EQUAL(   the_ider.get_id_of_name( string_ref{ emptiness_str  } ), 1                            );
	BOOST_CHECK      ( ! the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 2 );

	the_ider.clear();

	BOOST_CHECK      (   the_ider.empty()   );
	BOOST_CHECK_EQUAL(   the_ider.size(), 0 );
}


BOOST_AUTO_TEST_SUITE(largest_number_if_names_all_numeric_integers_fn)

BOOST_AUTO_TEST_CASE(returns_largest_negative) {
	id_of_str_bidirnl the_ider;

	the_ider.add_name( "-6"s );
	the_ider.add_name( "-7"s );
	the_ider.add_name( "-9"s );
	the_ider.add_name( "-8"s );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), make_optional( static_cast<ptrdiff_t>( -6 ) ) );
}

BOOST_AUTO_TEST_CASE(handles_non_numeric) {
	id_of_str_bidirnl the_ider;

	the_ider.add_name( "-6"s  );
	the_ider.add_name( "-7"s  );
	the_ider.add_name( "-9"s  );
	the_ider.add_name( "-8"s  );
	the_ider.add_name( "bob"s );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), none );
}


BOOST_AUTO_TEST_CASE(rejects_scientific_notation) {
	id_of_str_bidirnl the_ider;

	the_ider.add_name( "1e2"s  );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), none );
}

BOOST_AUTO_TEST_CASE(rejects_empty) {
	id_of_str_bidirnl the_ider;

	the_ider.add_name( ""s  );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), none );
}

BOOST_AUTO_TEST_CASE(rejects_single_dash) {
	id_of_str_bidirnl the_ider;

	the_ider.add_name( "-"s  );

	BOOST_CHECK_EQUAL( largest_number_if_names_all_numeric_integers( the_ider ), none );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()
