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

#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"
#include "options/options_block/data_dirs_options_block.h"
#include "structure/geometry/coord_list.h"
#include "superposition/superposition_context.h"
#include "test/global_test_constants.h"

using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::opts;
using namespace cath::sup;

using namespace std;

namespace cath {
	namespace test {

		/// \brief The superposition_context_test_suite_fixture to assist in testing superposition_context
		struct superposition_context_test_suite_fixture : protected global_test_constants {
		protected:
			~superposition_context_test_suite_fixture() noexcept = default;

			const coord_list            coord_list_1{ { coord{  1.0,  0.0,  0.0 }, coord{  2.0,   0.0,   0.0 } } };
			const coord_list            coord_list_2{ { coord{  0.0, -1.0,  0.0 }, coord{  0.0,  -2.0,   0.0 } } };
			const pdb_list              pdbs{ pdb_vec{ 2, pdb{} } };
			const str_vec               names{ "1c0pA01", "1hdoA00" };
			const superposition         the_sup{ create_pairwise_superposition( coord_list_1, coord_list_2 ) };
			const superposition_context the_sup_con{ pdbs, names, the_sup };
			const string                json_string = R"({"entries":[{"name":"1c0pA01","transformation":{"translation":)"
				R"({"x":"0","y":"0","z":"0"},)"
				R"("rotation":)"
				R"([["1","0","0"],)"
				R"(["0","1","0"],)"
				R"(["0","0","1"])"
				R"(]}},{"name":"1hdoA00","transformation":{"translation":)"
				R"({"x":"0","y":"0","z":"0"},)"
				R"("rotation":)"
				R"([["0","-1","0"],)"
				R"(["0","0","-1"],)"
				R"(["1","0","0"])"
				R"(]}}]})" "\n";
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(superposition_context_test_suite, cath::test::superposition_context_test_suite_fixture)

BOOST_AUTO_TEST_CASE(load_pdbs_from_names_copy_leaves_orig_empty_pdbs) {
	const auto loaded_sup_con = load_pdbs_from_names_copy( the_sup_con, build_data_dirs_options_block_of_dir( TEST_SOURCE_DATA_DIR() ) );
	BOOST_REQUIRE_EQUAL( the_sup_con.get_pdbs_cref().size(), 2 );
	BOOST_CHECK_EQUAL  ( the_sup_con.get_pdbs_cref()[ 0 ].get_num_residues(), 0 );
	BOOST_CHECK_EQUAL  ( the_sup_con.get_pdbs_cref()[ 1 ].get_num_residues(), 0 );
}

BOOST_AUTO_TEST_CASE(load_pdbs_from_names_copy_sets_two_pdbs_with_199_and_205_residues) {
	const auto loaded_sup_con = load_pdbs_from_names_copy( the_sup_con, build_data_dirs_options_block_of_dir( TEST_SOURCE_DATA_DIR() ) );
	BOOST_REQUIRE_EQUAL( loaded_sup_con.get_pdbs_cref().size(), 2 );
	BOOST_CHECK_EQUAL  ( loaded_sup_con.get_pdbs_cref()[ 0 ].get_num_residues(), 199 );
	BOOST_CHECK_EQUAL  ( loaded_sup_con.get_pdbs_cref()[ 1 ].get_num_residues(), 205 );
}

BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_example_sup_con) {
	BOOST_CHECK_EQUAL(
		to_json_string( the_sup_con, false ),
		json_string
	);
}

BOOST_AUTO_TEST_CASE(from_json_string_works) {
	const auto from_json_string = superposition_context_from_json_string( json_string );
	BOOST_REQUIRE_EQUAL     ( from_json_string.get_pdbs_cref().size(),   2       );
	BOOST_CHECK_EQUAL       ( from_json_string.get_pdbs_cref()[ 0 ].get_num_residues(), 0 );
	BOOST_CHECK_EQUAL       ( from_json_string.get_pdbs_cref()[ 1 ].get_num_residues(), 0 );

	BOOST_CHECK_EQUAL_RANGES( from_json_string.get_names_cref(),         names   );

	BOOST_CHECK_EQUAL       ( from_json_string.get_superposition_cref(), the_sup );

	BOOST_CHECK_EQUAL       ( from_json_string.has_alignment(),          false   );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

