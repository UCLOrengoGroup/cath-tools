/// \file
/// \brief The bioplib_pdb class

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

#include <boost/filesystem.hpp>

#include "common/boost_addenda/test/boost_check_no_throw_diag.hpp"
#include "test/global_test_constants.hpp"
#include "structure/bioplib_facade/bioplib_pdb.hpp"
#include "structure/geometry/coord.hpp"

#include <string>

using namespace boost::filesystem;
using namespace cath;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The bioplib_pdb_test_suite_fixture to assist in testing bioplib_pdb
		struct bioplib_pdb_test_suite_fixture : protected global_test_constants {
		protected:
			~bioplib_pdb_test_suite_fixture() noexcept = default;

		public:

			void check_chain_label(const bioplib_pdb &arg_bioplib_pdb,
			                       const char &arg_char
			                       ) {
				const PDB *  ptr    = arg_bioplib_pdb.get_ptr();
				const size_t natoms = arg_bioplib_pdb.get_natoms();
				for (size_t atom_ctr = 0; atom_ctr < natoms; ++atom_ctr) {
					BOOST_REQUIRE( ptr != nullptr );
					assert( ptr != nullptr ); // Post-BOOST_REQUIRE() assert() to appease clang's analyzer
					BOOST_CHECK_EQUAL( string{ arg_char }, string( ptr->chain ) );
					ptr = ptr->next;
				}
			}
			void check_numbering(const bioplib_pdb &arg_bioplib_pdb,
			                     const char &arg_char
			                     ) {
				const PDB *  ptr    = arg_bioplib_pdb.get_ptr();
				const size_t natoms = arg_bioplib_pdb.get_natoms();
				for (size_t atom_ctr = 0; atom_ctr < natoms; ++atom_ctr) {
					BOOST_REQUIRE( ptr != nullptr );
					assert( ptr != nullptr ); // Post-BOOST_REQUIRE() assert() to appease clang's analyzer
					BOOST_CHECK_EQUAL( string{ arg_char }, string(ptr->chain));
					ptr = ptr->next;
				}
			}
			void check_natoms(const bioplib_pdb &arg_bioplib_pdb,
			                  const size_t      &arg_natoms
			                  ) {
				BOOST_CHECK_EQUAL(arg_bioplib_pdb.get_natoms(), arg_natoms);
			}
		};
	}
}

BOOST_FIXTURE_TEST_SUITE(bioplib_pdb_test_suite, cath::test::bioplib_pdb_test_suite_fixture)

/// \brief Check that bioplib_pdb's set_chain_label() correctly sets the chain label
BOOST_AUTO_TEST_CASE(set_chain_label) {
	bioplib_pdb my_pdb;
	my_pdb.read_file(EXAMPLE_A_PDB_FILENAME().string());
	my_pdb.set_chain_label( chain_label( 'Z' ) );
	check_chain_label(my_pdb, 'Z');
	my_pdb.set_chain_label( chain_label( 'F') );
	check_chain_label(my_pdb, 'F');
}

/// \brief Check that bioplib_pdb's rotate() correctly rotates the coordinates
///
/// \todo Make this test check the result of rotating a bioplib_pdb
BOOST_AUTO_TEST_CASE(rotate) {
	bioplib_pdb my_pdb;
	my_pdb.read_file(EXAMPLE_A_PDB_FILENAME().string());
	BOOST_CHECK( true );
}

/// \brief Check that bioplib_pdb's operator+() correctly rotates the coordinates
///
/// \todo Make this test check the result of rotating a bioplib_pdb
BOOST_AUTO_TEST_CASE(translate) {
	bioplib_pdb my_pdb;
	my_pdb.read_file(EXAMPLE_A_PDB_FILENAME().string());
	const coord coord_123(1.0, 2.0, 3.0);
	my_pdb += coord_123;
	my_pdb += coord_123;
	my_pdb -= coord_123;
	BOOST_CHECK( true );
}

/// \brief Check that the copy ctor does not throw (and can copy-construct a const from a const)
BOOST_AUTO_TEST_CASE(copy_constructed_does_not_throw) {
	bioplib_pdb source_pdb;
	source_pdb.read_file(EXAMPLE_A_PDB_FILENAME().string());
	BOOST_REQUIRE_NO_THROW(const bioplib_pdb dest_pdb1(source_pdb));
	const bioplib_pdb dest_pdb1(source_pdb);
	BOOST_CHECK_NO_THROW_DIAG(const bioplib_pdb dest_pdb2(dest_pdb1));
}

/// \brief Check that the copy ctor produces a separate object
///
/// Running Valgrind on this testcase should catch any memory mischief
BOOST_AUTO_TEST_CASE(copy_constructed_is_independent) {
	bioplib_pdb source_pdb;
	source_pdb.read_file(EXAMPLE_A_PDB_FILENAME().string());
	bioplib_pdb dest_pdb(source_pdb);
	source_pdb.set_chain_label( chain_label( 'Z') );
	check_chain_label(source_pdb, 'Z');
	check_chain_label(dest_pdb,   'A');

	dest_pdb.set_chain_label( chain_label( 'F') );
	check_chain_label(source_pdb, 'Z');
	check_chain_label(dest_pdb,   'F');
}

BOOST_AUTO_TEST_SUITE_END()
