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

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "common/boost_check_no_throw_diag.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"
#include "test/global_test_constants.h"
#include "test/log_to_ostream_guard.h"

#include <vector>

using namespace boost::algorithm;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::algorithm::icontains;
using boost::algorithm::join;

namespace cath {
	namespace test {

		/// \brief The pdb_test_suite_fixture to assist in testing pdb
		struct pdb_test_suite_fixture : protected global_test_constants {
		protected:
			~pdb_test_suite_fixture() noexcept = default;

			/// \brief Check that the number of atoms in a vector of PDBs matches the expected numbers
			void check_nums_of_atoms(const pdb_list &,
			                         const size_vec &) const;
		};

	}
}

/// \brief Check that the number of atoms in a vector of PDBs matches the expected numbers
void cath::test::pdb_test_suite_fixture::check_nums_of_atoms(const pdb_list &arg_pdbs,                  ///< TODOCUMENT
                                                             const size_vec &arg_expected_nums_of_atoms ///< TODOCUMENT
                                                             ) const {
	// Grab the numbers of atoms
	size_vec got_nums_of_atoms;
	got_nums_of_atoms.reserve(arg_pdbs.size());
	for (const pdb &arg_pdb : arg_pdbs) {
		got_nums_of_atoms.push_back( arg_pdb.get_num_atoms() );
	}

	// Check that they match
	BOOST_CHECK_EQUAL_RANGES( arg_expected_nums_of_atoms, got_nums_of_atoms );
}

BOOST_FIXTURE_TEST_SUITE(pdb_test_suite, cath::test::pdb_test_suite_fixture)

/// \brief Test the parsing of END separated pdbs by sticking every combination of 0-3 ENDS before, between and after two 1-atom PDBs
///        (except there must be at least one END between)
BOOST_AUTO_TEST_CASE(check_parsing_of_end_separators) {
	const size_t MAX_NUM_CONSECUTIVE_ENDS(3);
	for (size_t num_before_ends = 0; num_before_ends < MAX_NUM_CONSECUTIVE_ENDS; ++num_before_ends) {
		for (size_t num_between_ends = 1; num_between_ends < MAX_NUM_CONSECUTIVE_ENDS; ++num_between_ends) {
			for (size_t num_after_ends = 0; num_after_ends < MAX_NUM_CONSECUTIVE_ENDS; ++num_after_ends) {
				stringstream end_separated_stream;
				end_separated_stream << join(str_vec(num_before_ends,  "END"), "\n") << endl;
				end_separated_stream << "ATOM   2952  OXT ALA   385      70.681 -13.748  36.367  1.00 26.84           O" << endl;
				end_separated_stream << join(str_vec(num_between_ends, "END"), "\n") << endl;
				end_separated_stream << "ATOM   2968  OXT ARG   387      20.593  77.271 -19.667  1.00  0.00           O" << endl;
				end_separated_stream << join(str_vec(num_after_ends,   "END"), "\n") << endl;
				check_nums_of_atoms(
					read_end_separated_pdb_files(end_separated_stream),
					{ 1_z, 1_z }
				);
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(parses_mse_resiude_as_unk_x) {
	// Example is from chain B of 4c9a
	const string input_string = R"(ATOM   1861  N   THR B 138     -33.417  42.721 103.639  1.00142.96           N  
ATOM   1862  CA  THR B 138     -33.726  41.382 103.094  1.00124.05           C  
ATOM   1863  C   THR B 138     -34.247  41.593 101.682  1.00122.85           C  
ATOM   1864  O   THR B 138     -34.642  42.709 101.356  1.00126.34           O  
ATOM   1865  CB  THR B 138     -34.816  40.679 103.927  1.00108.88           C  
ATOM   1866  N   MSE B 139     -34.243  40.568 100.830  1.00114.60           N  
ATOM   1867  CA  MSE B 139     -34.924  40.699  99.523  1.00114.24           C  
ATOM   1868  C   MSE B 139     -36.170  39.811  99.512  1.00107.49           C  
ATOM   1869  O   MSE B 139     -36.514  39.225  98.475  1.00116.82           O  
ATOM   1870  CB  MSE B 139     -33.991  40.377  98.339  1.00105.83           C  
ATOM   1871  N   VAL B 140     -36.867  39.765 100.654  1.00 94.64           N  
ATOM   1872  CA  VAL B 140     -37.904  38.745 100.917  1.00101.80           C  
ATOM   1873  C   VAL B 140     -39.124  39.279 101.672  1.00100.30           C  
ATOM   1874  O   VAL B 140     -39.041  40.306 102.330  1.00 97.38           O  
ATOM   1875  CB  VAL B 140     -37.331  37.579 101.757  1.00 93.10           C)";
	istringstream input_ss{ input_string };
	const pdb the_pdb = read_pdb_file( input_ss );
	BOOST_REQUIRE_EQUAL( the_pdb.get_num_residues(), 3 );
	BOOST_CHECK_EQUAL( get_amino_acid( the_pdb.get_residue_cref_of_index__backbone_unchecked( 1 ) ).get_code(), "UNK");
	// \todo: Create a better way to refer to check for UNK in code/tests without a string literal
}


BOOST_AUTO_TEST_CASE(skips_and_warns_on_silly_amino_acid) {
	// Example is a modified version of chain B of 4c9a
	const string input_string = R"(ATOM   1861  N   THR B 138     -33.417  42.721 103.639  1.00142.96           N  
ATOM   1862  CA  THR B 138     -33.726  41.382 103.094  1.00124.05           C  
ATOM   1863  C   THR B 138     -34.247  41.593 101.682  1.00122.85           C  
ATOM   1864  O   THR B 138     -34.642  42.709 101.356  1.00126.34           O  
ATOM   1865  CB  THR B 138     -34.816  40.679 103.927  1.00108.88           C  
ATOM   1866  N   FOO B 139     -34.243  40.568 100.830  1.00114.60           N  
ATOM   1867  CA  FOO B 139     -34.924  40.699  99.523  1.00114.24           C  
ATOM   1868  C   FOO B 139     -36.170  39.811  99.512  1.00107.49           C  
ATOM   1869  O   FOO B 139     -36.514  39.225  98.475  1.00116.82           O  
ATOM   1870  CB  FOO B 139     -33.991  40.377  98.339  1.00105.83           C  
ATOM   1871  N   VAL B 140     -36.867  39.765 100.654  1.00 94.64           N  
ATOM   1872  CA  VAL B 140     -37.904  38.745 100.917  1.00101.80           C  
ATOM   1873  C   VAL B 140     -39.124  39.279 101.672  1.00100.30           C  
ATOM   1874  O   VAL B 140     -39.041  40.306 102.330  1.00 97.38           O  
ATOM   1875  CB  VAL B 140     -37.331  37.579 101.757  1.00 93.10           C)";
	ostringstream test_ss;
	istringstream input_ss{ input_string };
	const log_to_ostream_guard the_guard{ test_ss };

	const pdb the_pdb = read_pdb_file( input_ss );

	BOOST_REQUIRE_EQUAL( the_pdb.get_num_residues(), 2 );
	BOOST_CHECK( icontains( test_ss.str(), "skip" ) );
}
BOOST_AUTO_TEST_SUITE_END()
