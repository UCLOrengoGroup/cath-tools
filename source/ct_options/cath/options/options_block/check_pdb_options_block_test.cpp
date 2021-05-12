/// \file
/// \brief The check_pdb_options_block test suite

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

#include <filesystem>

#include <boost/test/unit_test.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/file/temp_file.hpp"
#include "cath/options/options_block/check_pdb_options_block.hpp"
#include "cath/options/options_block/options_block_tester.hpp"

using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

using ::std::filesystem::path;

namespace cath {
	namespace test {

		/// \brief The check_pdb_options_block_test_suite_fixture to assist in testing check_pdb_options_block
		struct check_pdb_options_block_test_suite_fixture : protected options_block_tester {
		protected:
			~check_pdb_options_block_test_suite_fixture() noexcept = default;

			check_pdb_options_block the_options_block;
			const string IGNORE_OPT = { "positional-that-should-be-ignored" };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(check_pdb_options_block_test_suite, cath::test::check_pdb_options_block_test_suite_fixture)

/// \brief Check that if permit_atoms is requested, then that is handled correctly
BOOST_AUTO_TEST_CASE(handles_permit_atoms) {
	const temp_file temp_file("cath_tools_test_temp_file.check_pdb_options_block.%%%%");
	const path      temp_file_filename = get_filename( temp_file );
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT,
		  "--" + check_pdb_options_block::PO_PERMIT,
		  // "--" + check_pdb_options_block::PO_PDB_FILE,
		  temp_file_filename.string() }
	);
	BOOST_CHECK_EQUAL( path(), the_options_block.get_pdb_file()        );
	BOOST_CHECK_EQUAL( true,   the_options_block.get_permit_no_atoms() );
}

/// \brief Check that if an existent PDB file is requested then that is handled correctly
BOOST_AUTO_TEST_CASE(handles_existent_pdb_file) {
	const temp_file temp_file("cath_tools_test_temp_file.check_pdb_options_block.%%%%");
	const path      temp_file_filename = get_filename( temp_file );
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT,
		  "--" + check_pdb_options_block::PO_PDB_FILE,
		  temp_file_filename.string() }
	);
	BOOST_CHECK_EQUAL( temp_file_filename, the_options_block.get_pdb_file()        );
	BOOST_CHECK_EQUAL( false,              the_options_block.get_permit_no_atoms() );
}

/// \brief Check that if a non-existent PDB file is requested then that is handled correctly
BOOST_AUTO_TEST_CASE(handles_non_existent_pdb_file) {
	parse_into_options_block(
		the_options_block,
		{ IGNORE_OPT,
		  "--"+check_pdb_options_block::PO_PDB_FILE,
		  "/some/file/that/does/not/exist" }
	);
	BOOST_CHECK_EQUAL( path("/some/file/that/does/not/exist"), the_options_block.get_pdb_file()        );
	BOOST_CHECK_EQUAL( false,                                  the_options_block.get_permit_no_atoms() );
}

BOOST_AUTO_TEST_SUITE_END()

