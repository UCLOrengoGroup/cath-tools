/// \file


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
#include <boost/lexical_cast.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/output_test_stream.hpp>

#include "common/boost_check_no_throw_diag.h"
#include "common/file/open_fstream.h"
#include "common/size_t_literal.h"
#include "common/test_predicate/files_equal.h"
#include "exception/runtime_error_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "test/global_test_constants.h"
#include "structure/bioplib_facade/bioplib_pdb.h"
#include "structure/geometry/coord.h"

#include <fstream>
#include <iostream>

using namespace boost::filesystem;
using namespace boost::test_tools;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::test;
using namespace std;

using boost::lexical_cast;

namespace cath {
	namespace test {

		/// \brief The pdb_base_test_suite_fixture to assist in testing pdb_base
		struct pdb_base_test_suite_fixture : public global_test_constants {
		protected:
			/// \brief The directory containing all the PDBs
			const path ALL_PDBS_DIR{ path("/") / "cath" / "data" / "current" / "pdb" };

			pdb_base_test_suite_fixture();
			~pdb_base_test_suite_fixture() noexcept;

			/// \brief Specify that the copy-ctor shouldn't be used
			pdb_base_test_suite_fixture(const pdb_base_test_suite_fixture &) = delete;
			/// \brief Specify that the copy-assign shouldn't be used
			pdb_base_test_suite_fixture & operator=(const pdb_base_test_suite_fixture &) = delete;
		};

	}
}

/// \brief TODOCUMENT
cath::test::pdb_base_test_suite_fixture::pdb_base_test_suite_fixture() {
	BOOST_REQUIRE( ! exists( path( MODIFIED_PDB_FILENAME() ) ) );
}

/// \brief TODOCUMENT
cath::test::pdb_base_test_suite_fixture::~pdb_base_test_suite_fixture() noexcept {
	try {
		if ( exists( path( MODIFIED_PDB_FILENAME() ) ) ) {
			 remove( path( MODIFIED_PDB_FILENAME() ) );
		}
	}
	catch (...) {
	}
}

namespace cath {
	namespace test {

		/// \brief A functor for reading and writing
		///
		/// The operator() method is a template so that it can be called with an object of
		/// each of the types in all_pdb_types using boost::mpl::foreach<>().
		///
		/// This copies the source file to the TEST_OUTPUT_DIRECTORY so that when
		/// TEST_OUTPUT_DIRECTORY is set to /dev/shm, the file only need be read from disk once.
		class pdb_read_write_comparer {
		private:
			const path orig_source_pdb;
			const path test_dir;
			const path test_source_pdb;
			path_vec compare_files;

		public:
			pdb_read_write_comparer(const path &arg_source_pdb
			                        ) : orig_source_pdb(arg_source_pdb),
			                            test_dir(pdb_base_test_suite_fixture::TEST_OUTPUT_DIRECTORY()), // This is more robust
		//	                            test_dir( path() / "dev" / "shm" ), // This is a bit naughty but is very fast for full_read_write_comparison
			                            test_source_pdb(test_dir / orig_source_pdb.filename()) {
				if (exists(test_source_pdb)) {
					remove(test_source_pdb);
				}
				copy_file(orig_source_pdb, test_source_pdb);
			}
			virtual ~pdb_read_write_comparer() {
				try {
					for (const path &compare_file : compare_files) {
						if (exists(compare_file)) {
							remove(compare_file);
						}
					}
					if (exists(test_source_pdb)) {
						remove(test_source_pdb);
					}
				}
				catch (...) { /// Prevent the destructor throwing any exceptions
				}
			}
			template <typename P> void operator()(P pdb_obj) {
				pdb_obj.read_file(test_source_pdb.string());
				const string test_output_file_suffix = "_test_rw_" + lexical_cast<string>(compare_files.size() + 1);
				const path test_output_file = (test_dir / (path(test_source_pdb.filename()).string() + test_output_file_suffix));
				if (exists(test_output_file)) {
					remove(test_output_file);
				}
				pdb_obj.append_to_file(test_output_file.string());
				compare_files.push_back(test_output_file);
			}
			void compare_all_files() {
		//		bool any_files_differed = false;
				if ( ! compare_files.empty() ) {
					// Loop over each of the files after the first and compare them with the first
					for (size_t file_ctr = 1; file_ctr < compare_files.size(); ++file_ctr) {
						BOOST_CHECK_FILES_EQUAL( compare_files.front(), compare_files[ file_ctr ] );
					}
				}
		//		return !any_files_differed;
			}
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(pdb_base_test_suite, pdb_base_test_suite_fixture)

/// \brief A type-list containing pairs of matching comparator functor types.
///
/// Note that I attempted to implement this without type-lists by storing
/// binary_function<> objects but this was because I had forgotten
/// the lesson that binary_function<> doesn't actually contain a declaration for
/// the operator() method (pure virtual or otherwise).
///
/// This means that you can't call operator() on a binary_function reference.
using all_pdb_types = boost::mpl::list<pdb, bioplib_pdb>;

/// \brief Check that the constructor and destructor of bioplib_pdb don't throw
BOOST_AUTO_TEST_CASE_TEMPLATE(ctor_and_dtor_does_not_throw, pdb_type, all_pdb_types) {
	BOOST_CHECK_NO_THROW_DIAG(pdb_type my_pdb);
}

/// \brief Check that bioplib_pdb's read_file() throws when given a non-existent file
BOOST_AUTO_TEST_CASE_TEMPLATE(read_non_existent_file_throws, pdb_type, all_pdb_types) {
	pdb_type my_pdb;
	BOOST_CHECK_THROW(my_pdb.read_file(NONEXISTENT_FILE()), runtime_error_exception);
}

/// \brief Check that bioplib_pdb's read_file() doesn't throw (for a sensible input file)
BOOST_AUTO_TEST_CASE_TEMPLATE(read_file_does_not_throw, pdb_type, all_pdb_types) {
	pdb_type my_pdb;
	BOOST_CHECK_NO_THROW_DIAG(my_pdb.read_file(EXAMPLE_A_PDB_FILENAME()));
	BOOST_CHECK_EQUAL(1532_z, my_pdb.get_num_atoms());
}

/// \brief Check that bioplib_pdb's read_file() doesn't throw (for a sensible input file)
BOOST_AUTO_TEST_CASE_TEMPLATE(read_modify_write_pdb, pdb_type, all_pdb_types) {
	pdb_type my_pdb;
	BOOST_CHECK_NO_THROW_DIAG(my_pdb.read_file(EXAMPLE_A_PDB_FILENAME()));
	my_pdb += coord(1.020394867, 2.0230697, 3.0239576);
	my_pdb.append_to_file(TEST_OUTPUT_FILENAME());
}

/// \brief Check that all_pdb_types produce the same result if asked to read and write an example PDB file
BOOST_AUTO_TEST_CASE(read_write_example_pdb) {
	pdb_read_write_comparer the_comparer(EXAMPLE_A_PDB_FILENAME());
	boost::mpl::for_each<all_pdb_types>( ref( the_comparer ) );
	the_comparer.compare_all_files();
}

/// \brief A VERY LONG test to compare reading and writing ALL PDBs
//BOOST_AUTO_TEST_CASE(full_read_write_comparison) {
//	const path all_pdbs_dir_path(ALL_PDBS_DIR);
//	directory_iterator my_dir_itr(all_pdbs_dir_path);
//	const directory_iterator end_of_dir;
//	for (const path &pdb_file : boost::iterator_range<directory_iterator>( my_dir_itr, end_of_dir ) ) {
//		if (path(pdb_file.filename()).string().length() == 4) {
//			cerr << "TEMPORARILY SKIPPING " << pdb_file << endl;
//			continue;
//		}
//		cerr << "Testing " << pdb_file << endl;
//		pdb_read_write_comparer the_comparer(pdb_file);
//		boost::mpl::for_each<all_pdb_types>(ref(the_comparer));
//		the_comparer.compare_all_files();
//	}
//}

// grep -P --colour 'A\s+(22|25)\s+' /cath/data/current/pdb/1cbnA

// \grep conflict /tmp/build_test_cath_ssap.out -B1 | \grep Testing | sort | \grep -Po '\d\w{3}\S\d{0,2}' | xargs -iVAR wc -l /cath/data/current/pdb/VAR | sort -g
// 69 /cath/data/current/pdb/2z70B
// 89 /cath/data/current/pdb/3l8lB
// 89 /cath/data/current/pdb/3l8lD
// [...]

// Testing /cath/data/current/pdb/3l8lB
// WARNING: Amino acid "TRP" for residue "11" on chain 'B' conflicts with previous amino acid "PHE"

BOOST_AUTO_TEST_SUITE_END()
