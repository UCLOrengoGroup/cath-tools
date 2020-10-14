/// \file
/// \brief The get_determinant test suite

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

#include <boost/test/unit_test.hpp>

#include "cath/common/gsl/get_determinant.hpp"
#include "cath/common/gsl/gsl_matrix_wrp.hpp"

// using namespace cath::geom;
using namespace cath::geom::detail;

BOOST_AUTO_TEST_SUITE(get_determinant_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	gsl_matrix_wrp bob{ 3, 3 };

	gsl_matrix_set( bob.get_ptr(), 0, 0, -2.0 );
	gsl_matrix_set( bob.get_ptr(), 0, 1,  2.0 );
	gsl_matrix_set( bob.get_ptr(), 0, 2, -3.0 );

	gsl_matrix_set( bob.get_ptr(), 1, 0, -1.0 );
	gsl_matrix_set( bob.get_ptr(), 1, 1,  1.0 );
	gsl_matrix_set( bob.get_ptr(), 1, 2,  3.0 );

	gsl_matrix_set( bob.get_ptr(), 2, 0,  2.0 );
	gsl_matrix_set( bob.get_ptr(), 2, 1,  0.0 );
	gsl_matrix_set( bob.get_ptr(), 2, 2, -1.0 );

	BOOST_CHECK_EQUAL( get_determinant( bob ), 18 );
}

BOOST_AUTO_TEST_SUITE_END()
