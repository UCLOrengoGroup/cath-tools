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
#include "file/options/data_dirs_options_block.h"
#include "test/superposition_fixture.h"

using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The superposition_context_test_suite_fixture to assist in testing superposition_context
		struct superposition_context_test_suite_fixture : protected superposition_fixture {
		protected:
			~superposition_context_test_suite_fixture() noexcept = default;
		};
	}
}

BOOST_FIXTURE_TEST_SUITE(superposition_context_test_suite, cath::test::superposition_context_test_suite_fixture)

BOOST_AUTO_TEST_CASE(load_pdbs_from_names_copy_leaves_orig_empty_pdbs) {
	const auto loaded_sup_con = load_pdbs_from_names_copy( the_sup_con, build_data_dirs_spec_of_dir( TEST_SOURCE_DATA_DIR() ) );
	BOOST_REQUIRE_EQUAL( the_sup_con.get_pdbs_cref().size(), 2 );
	BOOST_CHECK_EQUAL  ( the_sup_con.get_pdbs_cref()[ 0 ].get_num_residues(), 0 );
	BOOST_CHECK_EQUAL  ( the_sup_con.get_pdbs_cref()[ 1 ].get_num_residues(), 0 );
}

BOOST_AUTO_TEST_CASE(load_pdbs_from_names_copy_sets_two_pdbs_with_199_and_205_residues) {
	const auto loaded_sup_con = load_pdbs_from_names_copy( the_sup_con, build_data_dirs_spec_of_dir( TEST_SOURCE_DATA_DIR() ) );
	BOOST_REQUIRE_EQUAL( loaded_sup_con.get_pdbs_cref().size(), 2 );
	BOOST_CHECK_EQUAL  ( loaded_sup_con.get_pdbs_cref()[ 0 ].get_num_residues(), 199 );
	BOOST_CHECK_EQUAL  ( loaded_sup_con.get_pdbs_cref()[ 1 ].get_num_residues(), 205 );
}

BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_example_sup_con) {
	BOOST_CHECK_EQUAL(
		to_json_string( the_sup_con, false ),
		sup_context_json_str
	);
}

BOOST_AUTO_TEST_CASE(from_json_string_works) {
	const auto from_json_string = superposition_context_from_json_string( sup_context_json_str );
	BOOST_REQUIRE_EQUAL     ( from_json_string.get_pdbs_cref().size(),   2       );
	BOOST_CHECK_EQUAL       ( from_json_string.get_pdbs_cref()[ 0 ].get_num_residues(), 0 );
	BOOST_CHECK_EQUAL       ( from_json_string.get_pdbs_cref()[ 1 ].get_num_residues(), 0 );

	BOOST_CHECK_EQUAL_RANGES( from_json_string.get_names_cref(),         names   );

	BOOST_CHECK_EQUAL       ( from_json_string.get_superposition_cref(), the_sup );

	BOOST_CHECK_EQUAL       ( from_json_string.has_alignment(),          false   );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

