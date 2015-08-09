/// \file
/// \brief The ssap_scores_file_alignment_acquirer test suite

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


#include <boost/filesystem/path.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "alignment/alignment.h"
#include "alignment/io/alignment_io.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "options/acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer.h"
#include "superposition/superpose_orderer.h"
#include "test/global_test_constants.h"
#include "test/log_to_ostream_guard.h"

#include <iostream> /// ***** TEMPORARY *****

using namespace cath::file;
using namespace cath::opts;
using namespace std; /// ***** TEMPORARY? *****

using boost::filesystem::path;
using cath::log_to_ostream_guard;

namespace cath {
	namespace test {

		/// \brief The ssap_scores_file_alignment_acquirer_test_suite_fixture to assist in testing ssap_scores_file_alignment_acquirer
		struct ssap_scores_file_alignment_acquirer_test_suite_fixture : protected global_test_constants {
		protected:
			~ssap_scores_file_alignment_acquirer_test_suite_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(ssap_scores_file_alignment_acquirer_test_suite, cath::test::ssap_scores_file_alignment_acquirer_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	const path gluing_dir = "build-test-data/ssap_scores_alignment_gluing";
	const ssap_scores_file_alignment_acquirer the_acquirer( gluing_dir / "1o7iB00_3uljB00_4gs3A00.ssap_scores" );
	const pdb_list the_pdbs = read_pdb_files( {
		gluing_dir / "1o7iB00",
		gluing_dir / "3uljB00",
		gluing_dir / "4gs3A00",
	} );

	// Use a lambda to enclose the log_to_ostream_guard whilst allowing the result to be passed back out
	const auto aln_and_sptree = [&] () {
		ostringstream parse_ss;
		const log_to_ostream_guard parse_log_guard{ parse_ss };
		return the_acquirer.get_alignment_and_spanning_tree( the_pdbs );
	} ();

	BOOST_CHECK_EQUAL(
		alignment_as_fasta_string( aln_and_sptree.first, the_pdbs, { "1o7iB00", "3uljB00", "4gs3A00" } ),
		R"(>1o7iB00
-MEEKVGNLKPNMESVNVTVRVLEASEARQIQTKNGVRTISEAIVGDE--T-----GRVKLTLWGKHAG----SIKEGQVVKIENAWTTAF--KGQVQLNAGSKTKIAEASE------DGFPESSQIPENTPTA
>3uljB00
G----------DPQVLRG-GHCKWFNVM----------GFG--FISMTEGSPL--ENPVDVFVHQSKLYMGFRSLKEGEPVEFF-KS--S------KGFES----RVTGPGGNPCLGE----------------
>4gs3A00
-----------ENNTVTLVGKVFTPLEFSHEL---YGEKFFNFILEVP--RLSETKDYLPITISNRLFEG--MNLEVGTRVKIE-GQLRSYNRKLILTVFA---RDISVVPE----------------------
)"
		);
}
BOOST_AUTO_TEST_SUITE_END()

