/// \file
/// \brief The superposition_context test suite

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

#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"
#include "structure/geometry/coord_list.h"
#include "superposition/superposition_context.h"

using namespace cath::file;
using namespace cath::geom;
using namespace cath::sup;

namespace cath {
	namespace test {

		/// \brief The superposition_context_test_suite_fixture to assist in testing superposition_context
		struct superposition_context_test_suite_fixture {
		protected:
			~superposition_context_test_suite_fixture() noexcept = default;

			coord_list coord_list_1{ { coord{  1.0,  0.0,  0.0 }, coord{  2.0,   0.0,   0.0 } } };
			coord_list coord_list_2{ { coord{  0.0, -1.0,  0.0 }, coord{  0.0,  -2.0,   0.0 } } };
			superposition the_sup = create_pairwise_superposition( coord_list_1, coord_list_2 );
			superposition_context the_sup_con{
				pdb_list{ pdb_vec{ 2, pdb{} } },
				str_vec { "name_a", "name_b" },
				the_sup
			};
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(superposition_context_test_suite, cath::test::superposition_context_test_suite_fixture)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_example_sup_con) {
	BOOST_CHECK_EQUAL(
		to_json_string( the_sup_con, false ),
		R"({"entries":[{"name":"name_a","transformation":{"translation":)"
		R"({"x":"0","y":"0","z":"0"},)"
		R"("rotation":)"
		R"([["1","0","0"],)"
		R"(["0","1","0"],)"
		R"(["0","0","1"])"
		R"(]}},{"name":"name_b","transformation":{"translation":)"
		R"({"x":"0","y":"0","z":"0"},)"
		R"("rotation":)"
		R"([["0","-1","0"],)"
		R"(["0","0","-1"],)"
		R"(["1","0","0"])"
		R"(]}}]})" "\n"
	);
}

BOOST_AUTO_TEST_SUITE_END()

