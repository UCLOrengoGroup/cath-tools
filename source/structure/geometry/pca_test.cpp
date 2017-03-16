/// \file
/// \brief The pca test suite

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

#include "structure/geometry/coord.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/geometry/line.hpp"
#include "structure/geometry/pca.hpp"
#include "structure/structure_type_aliases.hpp"

using namespace cath::geom;

BOOST_AUTO_TEST_SUITE(pca_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	const coord_list example_coords{ coord_vec{
		coord{ 10.2, 2.0, -9.7 },
		coord{  5.2, 2.0, -4.7 },
		coord{ -1.8, 2.0,  2.3 },
		coord{ -2.8, 2.0,  3.3 },
		coord{ -9.8, 2.0, 10.3 },
	} };

	const auto lobf = line_of_best_fit( example_coords );

	BOOST_CHECK_EQUAL( lobf.get_point_on_line(), coord(  0.2,      2.0, 0.3      ) );
	BOOST_CHECK_EQUAL( lobf.get_dirn         (), coord( -0.707107, 0.0, 0.707107 ) );
}

BOOST_AUTO_TEST_SUITE_END()
