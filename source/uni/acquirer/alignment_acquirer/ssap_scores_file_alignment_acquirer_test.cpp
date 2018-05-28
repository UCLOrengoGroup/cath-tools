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

#include <boost/test/auto_unit_test.hpp>

#include "acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer.hpp"
#include "alignment/alignment.hpp"
#include "alignment/io/alignment_io.hpp"
#include "common/boost_addenda/log/stringstream_log_sink.hpp"
#include "file/strucs_context.hpp"
#include "test/global_test_constants.hpp"

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::align;
using namespace cath::file;
using namespace cath::test;

using std::ostringstream;

namespace cath {
	namespace test {

		/// \brief The ssap_scores_file_alignment_acquirer_test_suite_fixture to assist in testing ssap_scores_file_alignment_acquirer
		struct ssap_scores_file_alignment_acquirer_test_suite_fixture : protected global_test_constants {
		protected:
			~ssap_scores_file_alignment_acquirer_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(ssap_scores_file_alignment_acquirer_test_suite, ssap_scores_file_alignment_acquirer_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	const ssap_scores_file_alignment_acquirer the_acquirer( TEST_SSAP_ALIGNMENT_GLUING_DATA_DIR() / "1o7iB00_3uljB00_4gs3A00.ssap_scores" );
	const pdb_list the_pdbs = read_pdb_files( {
		TEST_SSAP_ALIGNMENT_GLUING_DATA_DIR() / "1o7iB00",
		TEST_SSAP_ALIGNMENT_GLUING_DATA_DIR() / "3uljB00",
		TEST_SSAP_ALIGNMENT_GLUING_DATA_DIR() / "4gs3A00",
	} );

	// Use a lambda to enclose the stringstream_log_sink whilst allowing the result to be passed back out
	const auto aln_and_sptree = [&] {
		const stringstream_log_sink log_sink;
		return the_acquirer.get_alignment_and_spanning_tree( strucs_context{ the_pdbs } );
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

// Alignment from fixing the SSAP scores so that gluing happens via 3uljB00:
// >1o7iB00
// MEEKVGNLKPNMESVNVT-VRVLEASEARQIQTKNGVR-------TIS--EAIVG---DE--TGRVKLTLWGKH----AGSIKEGQVVKIENAWTTAF--K------GQVQLNAGSKT-KIAEASEDGFPESSQIPENTPTA
// >3uljB00
// ---------GDPQVLRG--GHCKWFNVM-----------------GFG--FISMTEGSPL--ENPVDVFVHQSKLYMGFRSLKEGEPVEF--FK---S--S------K--GFES-----RVTGPGGNPCLGE----------
// >4gs3A00
// ----------ENNTVTL-VGKVFTPLEF----------SHELYGEKFFNFILEVP--RLSETKDYLPITISNRLFEG--MNLEVGTRVKI--EG---QLRSYNRKLIL--TVFA----RDISVVPE----------------

// Fully refined alignment:
// >1o7iB00
// MEEKVGNLKPNME-SVNVTVRVLE-ASEARQIQTKNGVRTISEAIVGDET-----GRVKLTLWGKHAG----S-IKEGQVVKIENAWTTAF--KGQVQLNAGSKTKIAEASEDGFPESSQIPENTPTA-
// >3uljB00
// ---------GSDPQVLRGSGHCKW-FNVRM---------GFGFISMTSREGSPLENPVDVFVHQSKLYMEGFRSLKEGEPVEFT-FKKSSK------GFES---LRVTGPGGNP----------CLGNE
// >4gs3A00
// ----------ENN-TVTLVGKVFTPLEFSHELY----GEKFFNFILEVPRLSETKDYLPITISNRLFEG---MNLEVGTRVKIE-GQLRSYNRKLILTVFA---RDISVVPE-----------------


// aln_glue_style::INCREMENTALLY_WITH_PAIR_REFINING
// >1o7iB00
// MEEKVGNLKPNMESVNVTVRVLE-ASEARQIQTKNGVRTISEAIVGDET-------GRVKLTLWGKHAG----S-IKEGQVVKIENAWTTAF--KGQVQLNAGSKTKIAEASEDGFPESSQIPENTPTA
// >3uljB00
// ------GSD--PQVLRGSGHCKW-FNVRM---------GFGF--ISMTSREGSPLENPVDVFVHQSKLYMEGFRSLKEGEPVEFT-FKKS-S-----KGFES---LRVTGPGGNPCLGNE---------
// >4gs3A00
// -------EN---NTVTLVGKVFTPLE-FSHE-L-Y-GEKFFNFILEVPR--LSETKDYLPITISNRLFEG---MNLEVGTRVKIE-GQLRSYNRKLILTVFA---RDISVVPE----------------


// Got :
// >1o7iB00
// MEEKVGNLKPNMESVNVTVRVLE-ASEARQIQTKNGVRTISEAIVGDE------T-GRVKLTLWGKHAG----S-IKEGQVVKIENAWTTAF--KGQVQLNAGSKTKIAEASEDGFPESSQIPENTPTA
// >3uljB00
// GS------D--PQVLRGSGHCKW-FNVRM---------GFGF--ISMTSREGSPLENPVDVFVHQSKLYMEGFRSLKEGEPVEFT-FKKSSK------GFES---LRVTGPGGNPCLGNE---------
// >4gs3A00
// -------EN--N-TVTLVGKVFTPLEFSHELY----GEKFFNFILEVP--RLSETKDYLPITISNRLFEG---MNLEVGTRVKIE-GQLRSYNRKLILTVFA---RDISVVPE----------------


		);
}

BOOST_AUTO_TEST_CASE(accepts_empty_scores_and_single_structure) {
	const ssap_scores_file_alignment_acquirer the_acquirer( TEST_SSAP_ALIGNMENT_GLUING_DATA_DIR() / "empty.ssap_scores" );
	const pdb_list the_pdbs = read_pdb_files( {
		TEST_SSAP_ALIGNMENT_GLUING_DATA_DIR() / "1o7iB00"
	} );

	// Use a lambda to enclose the stringstream_log_sink whilst allowing the result to be passed back out
	const auto aln_and_sptree = [&] {
		const stringstream_log_sink log_sink;
		return the_acquirer.get_alignment_and_spanning_tree( strucs_context{ the_pdbs } );
	} ();

	BOOST_CHECK_EQUAL(
		alignment_as_fasta_string( aln_and_sptree.first, the_pdbs, { "1o7iB00" } ),
		R"(>1o7iB00
MEEKVGNLKPNMESVNVTVRVLEASEARQIQTKNGVRTISEAIVGDETGRVKLTLWGKHAGSIKEGQVVKIENAWTTAFKGQVQLNAGSKTKIAEASEDGFPESSQIPENTPTA
)"
		);
}


BOOST_AUTO_TEST_SUITE_END()

