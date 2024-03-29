/// \file
/// \brief The residue_name test suite

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

#include "cath/biocore/residue_name.hpp"

using namespace ::cath;

BOOST_AUTO_TEST_SUITE(residue_name_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( to_string( residue_name(         ) ), "null_res" );
	BOOST_CHECK_EQUAL( to_string( residue_name( -5      ) ), "-5"       );
	BOOST_CHECK_EQUAL( to_string( residue_name( -5, 'A' ) ), "-5A"      );
}

BOOST_AUTO_TEST_CASE(negative_number_check_works) {
	static_assert(   has_strictly_negative_residue_number( residue_name( -1      ) ) );
	static_assert(   has_strictly_negative_residue_number( residue_name( -1, 'A' ) ) );

	static_assert( ! has_strictly_negative_residue_number( residue_name(  0      ) ) );
	static_assert( ! has_strictly_negative_residue_number( residue_name(  0, 'A' ) ) );

	static_assert( ! has_strictly_negative_residue_number( residue_name(  1      ) ) );
	static_assert( ! has_strictly_negative_residue_number( residue_name(  1, 'A' ) ) );

	static_assert( ! has_strictly_negative_residue_number( residue_name(         ) ) );

	BOOST_TEST( true );
}

BOOST_AUTO_TEST_SUITE_END()
