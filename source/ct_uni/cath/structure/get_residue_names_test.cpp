/// \file
/// \brief A test suite to check that getting a list of residue names from a pdb or a dssp works as expected

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

#include <boost/filesystem.hpp>

#include "cath/common/file/open_fstream.hpp"
#include "cath/file/dssp_wolf/dssp_file.hpp"
#include "cath/file/dssp_wolf/dssp_file_io.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"
#include "cath/test/global_test_constants.hpp"

#include <fstream>
#include <string>

namespace cath { namespace test { } }

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::test;
using namespace ::std;

using ::boost::filesystem::path;

/// \file
///
/// Commands used to produce the files:
///
///     cat 1x0pJ | perl -ni -e 'if ($_ =~ /^ATOM/) {print substr($_, 22, 5)."\n"; } ' | grep -Po '\S+' | uniq > 1x0pJ.residue_names
///     dssp-2.0.4-linux-amd64 1x0pJ > 1x0pJ.dssp
///     cat 1x0pJ.dssp | grep -PA99999 '^\s*#' | tail -n +2 | perl -ni -e 'print substr($_, 6, 5)."\n";' | grep -Po '\S+' | uniq > 1x0pJ.dssp_residue_names
///
/// Tests:
///  - 1bmv2: residues are out of order
///  - 1t63A: last residue has a high number and an insert code
///  - 1x0pJ: last residue has a high number

namespace cath {
	namespace test {

		/// \brief A simple enum for specifying whether the PDB code or the DSSP code is being tested
		enum class get_residue_ids_test_filetype : bool {
			PDB,
			DSSP
		};


		/// \brief The get_residue_names_test_suite_fixture to assist in testing getting residue names
		struct get_residue_names_test_suite_fixture : protected global_test_constants {
		private:
			const string pdb_extension                   = "";
			const string dssp_extension                  = ".dssp";
			const string correct_residue_names_extension = ".residue_names";

		protected:
			~get_residue_names_test_suite_fixture() noexcept = default;

		public:

			void check_get_residue_ids(const string                        &prm_chain_id,
			                           const get_residue_ids_test_filetype &prm_filetype
			                           ) const {
				residue_id_vec got_residue_ids;
				switch ( prm_filetype ) {
					case ( get_residue_ids_test_filetype::PDB ) : {
						const path pdb_file ( TEST_RESIDUE_IDS_DATA_DIR() / ( prm_chain_id + pdb_extension ) );
						pdb my_pdb;
						my_pdb.read_file( pdb_file.string() );
						got_residue_ids = my_pdb.get_residue_ids_of_first_chain__backbone_unchecked();
						break;
					}
					case ( get_residue_ids_test_filetype::DSSP ) : {
						const path the_dssp_file( TEST_RESIDUE_IDS_DATA_DIR() / ( prm_chain_id + dssp_extension ) );
						const dssp_file my_dssp = read_dssp_file( the_dssp_file );
						got_residue_ids = get_residue_ids( my_dssp, true );
						break;
					}
				}

				// Test the residues names against a file containing the correct answer
				//
				// Note: Not using output_test_stream to check the pattern file because that wouldn't complain if "got"
				//       was shorter than "expected".
				//
				// \todo I think I should write some better file/string matching tools (which could use some of the code from below):
				//  * default check that "got" length matches "expected" length
				//  * offer a bunch of check tools (eg equals, contains, regexp)
				//  * optionally load a pattern file as a vector of strings, one per line
				//  * provide better context information on mismatch
				const path correct_residue_names_file(TEST_RESIDUE_IDS_DATA_DIR() / (prm_chain_id + correct_residue_names_extension));
				ifstream expected_stream;
				open_ifstream( expected_stream, correct_residue_names_file );
				residue_id_vec expected_values;
				copy(
					istream_iterator<residue_id>( expected_stream ),
					istream_iterator<residue_id>(),
					back_inserter( expected_values )
				);
				expected_stream.close();

				BOOST_CHECK_EQUAL_RANGES( expected_values, got_residue_ids );
			}
		};

	}  // namespace test
}  // namespace cath

/// \brief A test suite to check that getting a list of residue names from a pdb or a dssp works as expected
///
/// ATM, there are test cases hard coding the six combinations of chain ID and filetype as a quick way to give them helpful names
/// \todo Improve this (read Boost test documentation (eg understand stuff in "Manually registered nullary function based test case")
////                    for the best way to do this)
BOOST_FIXTURE_TEST_SUITE(get_residue_names_test_suite, get_residue_names_test_suite_fixture)

/// \brief Test the reading of residue names from the pdb for 1bmv2
BOOST_AUTO_TEST_CASE(pdb_1bmv2) {
	check_get_residue_ids("1bmv2", get_residue_ids_test_filetype::PDB);
}

/// \brief Test the reading of residue names from the dssp for 1bmv2
BOOST_AUTO_TEST_CASE(dssp_1bmv2) {
	check_get_residue_ids("1bmv2", get_residue_ids_test_filetype::DSSP);
}

/// \brief Test the reading of residue names from the pdb for 1t63A
BOOST_AUTO_TEST_CASE(pdb_1t63A) {
	check_get_residue_ids("1t63A", get_residue_ids_test_filetype::PDB);
}

/// \brief Test the reading of residue names from the dssp for 1t63A
BOOST_AUTO_TEST_CASE(dssp_1t63A) {
	check_get_residue_ids("1t63A", get_residue_ids_test_filetype::DSSP);
}

/// \brief Test the reading of residue names from the pdb for 1x0pJ
BOOST_AUTO_TEST_CASE(pdb_1x0pJ) {
	check_get_residue_ids("1x0pJ", get_residue_ids_test_filetype::PDB);
}

/// \brief Test the reading of residue names from the dssp for 1x0pJ
BOOST_AUTO_TEST_CASE(dssp_1x0pJ) {
	check_get_residue_ids("1x0pJ", get_residue_ids_test_filetype::DSSP);
}

BOOST_AUTO_TEST_SUITE_END()
