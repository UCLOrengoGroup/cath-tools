///// \file
///// \brief The superposition I/O test suite

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

#include "superposition/superposition.h"
#include "superposition/io/superposition_io.h"

using namespace cath::geom;
using namespace cath::sup;

/// \brief The superposition_io_test_suite_fixture to assist in testing superposition I/O
struct superposition_io_test_suite_fixture {
protected:
	~superposition_io_test_suite_fixture() noexcept = default;

	coord_list coord_list_1{ { coord{  1.0,  0.0,  0.0 }, coord{  2.0,   0.0,   0.0 } } };
	coord_list coord_list_2{ { coord{  0.0, -1.0,  0.0 }, coord{  0.0,  -2.0,   0.0 } } };
	superposition the_sup = create_pairwise_superposition( coord_list_1, coord_list_2 );
};

BOOST_FIXTURE_TEST_SUITE(superposition_io_test_suite, superposition_io_test_suite_fixture)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_example_sup) {
	BOOST_CHECK_EQUAL(
		to_json_string( create_pairwise_superposition( coord_list_1, coord_list_2 ), false ),
		R"({"transformations":[{"translation":)"
		R"({"x":"0","y":"0","z":"0"},)"
		R"("rotation":)"
		R"([["1","0","0"],)"
		R"(["0","1","0"],)"
		R"(["0","0","1"])"
		R"(]},{"translation":)"
		R"({"x":"0","y":"0","z":"0"},)"
		R"("rotation":)"
		R"([["0","-1","0"],)"
		R"(["0","0","-1"],)"
		R"(["1","0","0"])"
		R"(]}]})"
		"\n"
	);
}

BOOST_AUTO_TEST_SUITE_END()

