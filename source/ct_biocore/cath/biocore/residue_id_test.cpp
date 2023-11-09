/// \file
/// \brief The residue_id test suite

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

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/biocore/residue_id.hpp"

using namespace ::cath;

BOOST_AUTO_TEST_SUITE(residue_id_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( to_string( make_residue_id( 'A'          ) ), "A:null_res" );
	BOOST_CHECK_EQUAL( to_string( make_residue_id( 'A', -5      ) ), "A:-5"       );
	BOOST_CHECK_EQUAL( to_string( make_residue_id( 'A', -5, 'A' ) ), "A:-5A"      );

	static_assert( make_residue_id( 'A'          ) == make_residue_id( 'A'          ) );
	static_assert( make_residue_id( 'A', -5      ) == make_residue_id( 'A', -5      ) );
	static_assert( make_residue_id( 'A', -5, 'A' ) == make_residue_id( 'A', -5, 'A' ) );
}

BOOST_AUTO_TEST_CASE(negative_number_check_works) {
	static_assert(   has_strictly_negative_residue_number( make_residue_id( 'A', -1      ) ) );
	static_assert(   has_strictly_negative_residue_number( make_residue_id( 'A', -1, 'A' ) ) );

	static_assert( ! has_strictly_negative_residue_number( make_residue_id( 'A',  0      ) ) );
	static_assert( ! has_strictly_negative_residue_number( make_residue_id( 'A',  0, 'A' ) ) );

	static_assert( ! has_strictly_negative_residue_number( make_residue_id( 'A',  1      ) ) );
	static_assert( ! has_strictly_negative_residue_number( make_residue_id( 'A',  1, 'A' ) ) );

	static_assert( ! has_strictly_negative_residue_number( make_residue_id( 'A'          ) ) );

	BOOST_TEST( true );
}

BOOST_AUTO_TEST_CASE(any_negative_number_check_works) {
	const residue_id_vec none_strictly_negative = { {
		make_residue_id( 'A',  0      ),
		make_residue_id( 'A',  0, 'A' ),
		make_residue_id( 'A',  1      ),
		make_residue_id( 'A',  1, 'A' ),
		make_residue_id( 'A'          ),
	} };

	BOOST_CHECK( ! has_any_strictly_negative_residue_numbers ( none_strictly_negative ) );

	const residue_id_vec one_strictly_negative = { {
		make_residue_id( 'A',  0      ),
		make_residue_id( 'A',  0, 'A' ),
		make_residue_id( 'A',  1      ),
		make_residue_id( 'A', -1, 'A' ),
		make_residue_id( 'A',  1, 'A' ),
		make_residue_id( 'A'          ),
	} };

	BOOST_CHECK(   has_any_strictly_negative_residue_numbers (  one_strictly_negative ) );
}

BOOST_AUTO_TEST_SUITE_END()
