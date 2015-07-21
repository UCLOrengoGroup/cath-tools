/// \file


/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "common/boost_check_no_throw_diag.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"
#include "test/global_test_constants.h"

#include <vector>

using namespace boost::algorithm;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::algorithm::join;

namespace cath {
	namespace test {

		/// \brief The pdb_test_suite_fixture to assist in testing pdb
		struct pdb_test_suite_fixture : public global_test_constants {
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

BOOST_AUTO_TEST_SUITE_END()
