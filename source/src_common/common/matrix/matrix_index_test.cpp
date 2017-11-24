/// \matrix
/// \brief The matrix_index test suite

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

#include <boost/test/auto_unit_test.hpp>

#include "common/matrix/matrix_index.hpp"

using namespace cath::common;

BOOST_AUTO_TEST_SUITE(matrix_index_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert( get_zero_index_of_strict_upper_half_matrix( 0, 1, 5 ) == 0, "" );
	static_assert( get_zero_index_of_strict_upper_half_matrix( 0, 2, 5 ) == 1, "" );
	static_assert( get_zero_index_of_strict_upper_half_matrix( 0, 3, 5 ) == 2, "" );
	static_assert( get_zero_index_of_strict_upper_half_matrix( 0, 4, 5 ) == 3, "" );

	static_assert( get_zero_index_of_strict_upper_half_matrix( 1, 2, 5 ) == 4, "" );
	static_assert( get_zero_index_of_strict_upper_half_matrix( 1, 3, 5 ) == 5, "" );
	static_assert( get_zero_index_of_strict_upper_half_matrix( 1, 4, 5 ) == 6, "" );

	static_assert( get_zero_index_of_strict_upper_half_matrix( 2, 3, 5 ) == 7, "" );
	static_assert( get_zero_index_of_strict_upper_half_matrix( 2, 4, 5 ) == 8, "" );

	static_assert( get_zero_index_of_strict_upper_half_matrix( 3, 4, 5 ) == 9, "" );

	BOOST_TEST( true );
}

BOOST_AUTO_TEST_SUITE_END()
