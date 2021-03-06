/// \file
/// \brief The extract_pdb_options_block test suite

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

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/file/temp_file.hpp"
#include "cath/options/options_block/extract_pdb_options_block.hpp"
#include "cath/options/options_block/options_block_tester.hpp"

using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

namespace {

	/// \brief The extract_pdb_options_block_test_suite_fixture to assist in testing extract_pdb_options_block
	struct extract_pdb_options_block_test_suite_fixture : protected options_block_tester {
	protected:
		~extract_pdb_options_block_test_suite_fixture() noexcept = default;

		extract_pdb_options_block the_options_block;
		const string IGNORE_OPT = { "positional-that-should-be-ignored" };
	};

} // namespace

BOOST_FIXTURE_TEST_SUITE(extract_pdb_options_block_test_suite, extract_pdb_options_block_test_suite_fixture)

// /// \brief Check that if permit_atoms is requested, then that is handled correctly
// BOOST_AUTO_TEST_CASE(handles_permit_atoms) {
// 	const temp_file temp_file("cath_tools_test_temp_file.extract_pdb_options_block.%%%%");
// 	const path      temp_file_filename = get_filename( temp_file );
// 	parse_into_options_block(
// 		the_options_block,
// 		{ IGNORE_OPT,
// 		  "--" + extract_pdb_options_block::PO_PERMIT,
// 		  // "--" + extract_pdb_options_block::PO_PDB_FILE,
// 		  temp_file_filename.string() }
// 	);
// 	BOOST_EXTRACT_EQUAL( path(), the_options_block.get_pdb_file()        );
// 	BOOST_EXTRACT_EQUAL( true,   the_options_block.get_permit_no_atoms() );
// }

// /// \brief Check that if an existent PDB file is requested then that is handled correctly
// BOOST_AUTO_TEST_CASE(handles_existent_pdb_file) {
// 	const temp_file temp_file("cath_tools_test_temp_file.extract_pdb_options_block.%%%%");
// 	const path      temp_file_filename = get_filename( temp_file );
// 	parse_into_options_block(
// 		the_options_block,
// 		{ IGNORE_OPT,
// 		  "--" + extract_pdb_options_block::PO_INPUT_PDB_FILE,
// 		  temp_file_filename.string() }
// 	);
// 	BOOST_EXTRACT_EQUAL( temp_file_filename, the_options_block.get_pdb_file()        );
// 	BOOST_EXTRACT_EQUAL( false,              the_options_block.get_permit_no_atoms() );
// }

// /// \brief Check that if a non-existent PDB file is requested then that is handled correctly
// BOOST_AUTO_TEST_CASE(handles_non_existent_pdb_file) {
// 	parse_into_options_block(
// 		the_options_block,
// 		{ IGNORE_OPT,
// 		  "--"+extract_pdb_options_block::PO_INPUT_PDB_FILE,
// 		  "/some/file/that/does/not/exist" }
// 	);
// 	BOOST_EXTRACT_EQUAL( path("/some/file/that/does/not/exist"), the_options_block.get_pdb_file()        );
// 	BOOST_EXTRACT_EQUAL( false,                                  the_options_block.get_permit_no_atoms() );
// }

BOOST_AUTO_TEST_SUITE_END()

