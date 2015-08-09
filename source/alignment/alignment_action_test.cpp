/// \file
/// \brief The alignment_action test suite

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

#include "alignment/alignment_action.h"

#include "alignment/alignment.h"
#include "alignment/io/alignment_io.h"
#include "common/test_predicate/istream_and_file_equal.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_io.h"
#include "structure/protein/protein_list.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"
#include "test/global_test_constants.h"

using namespace boost::filesystem;
using namespace cath;
using namespace cath::align;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The alignment_action_test_suite_fixture to assist in testing alignment_action
		struct alignment_action_test_suite_fixture : public global_test_constants {
		private:
			ostringstream test_stderr;

		protected:
			~alignment_action_test_suite_fixture() noexcept = default;

		public:
			alignment_action_test_suite_fixture();

			const path         root_dir           = { TEST_MULTI_SSAP_SUPERPOSE_DIR() };
			const path         expected_aln_file  = { root_dir / "expected_glue_output.fa" };
			const protein      protein_1g5aA03    = { read_protein_from_dssp_and_pdb( root_dir / "1g5aA03.dssp", root_dir / "1g5aA03", true, "1g5aA03" ) };
			const protein      protein_1r7aA02    = { read_protein_from_dssp_and_pdb( root_dir / "1r7aA02.dssp", root_dir / "1r7aA02", true, "1r7aA02" ) };
			const protein      protein_1wzaA02    = { read_protein_from_dssp_and_pdb( root_dir / "1wzaA02.dssp", root_dir / "1wzaA02", true, "1wzaA02" ) };
			const protein      protein_1zjaA02    = { read_protein_from_dssp_and_pdb( root_dir / "1zjaA02.dssp", root_dir / "1zjaA02", true, "1zjaA02" ) };
			const protein_list all_proteins       = { make_protein_list( { protein_1g5aA03, protein_1r7aA02, protein_1wzaA02, protein_1zjaA02 } ) };
			const alignment    aln_1wzaA02_1zjaA02= { read_alignment_from_cath_ssap_legacy_format( root_dir / "1wzaA021zjaA02.list", protein_1wzaA02, protein_1zjaA02, test_stderr ) };
			const alignment    aln_1g5aA03_1zjaA02= { read_alignment_from_cath_ssap_legacy_format( root_dir / "1g5aA031zjaA02.list", protein_1g5aA03, protein_1zjaA02, test_stderr ) };
			const alignment    aln_1g5aA03_1r7aA02= { read_alignment_from_cath_ssap_legacy_format( root_dir / "1g5aA031r7aA02.list", protein_1g5aA03, protein_1r7aA02, test_stderr ) };
		};

		/// \brief TODOCUMENT
		alignment_action_test_suite_fixture::alignment_action_test_suite_fixture() {
			if ( ! test_stderr.str().empty() ) {
				cerr << endl;
				cerr << endl;
				cerr << test_stderr.str();
				cerr << endl;
				cerr << endl;
			}
			assert( test_stderr.str().empty() );
		}
	}
}

BOOST_FIXTURE_TEST_SUITE(alignment_action_test_suite, cath::test::alignment_action_test_suite_fixture)

/// \brief Test that the code to glue alignments together works as expected
BOOST_AUTO_TEST_CASE(separately_glue_1wzaA02_1zjaA02_1g5aA03_1r7aA02) {
	const alignment aln_1g5aA03_1zjaA02_1wzaA02         = glue_two_alignments( aln_1g5aA03_1zjaA02,         1,
	                                                                           aln_1wzaA02_1zjaA02,         1  );
	const alignment aln_1g5aA03_1zjaA02_1wzaA02_1r7aA02 = glue_two_alignments( aln_1g5aA03_1zjaA02_1wzaA02, 0,
	                                                                           aln_1g5aA03_1r7aA02,         0  );

	// Check that the alignment matches the expected
	stringstream got_ss;
	write_alignment_as_fasta_alignment(
		got_ss,
		make_permuted_alignment( aln_1g5aA03_1zjaA02_1wzaA02_1r7aA02, { 0, 3, 2, 1 } ),
		all_proteins
	);
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( got_ss, "got_ss", expected_aln_file );
}

/// \brief Test that the code to automatically build an alignment from a spanning-tree works as expected
BOOST_AUTO_TEST_CASE(build_alignment_from_parts_works) {
	const alignment new_alignment = build_alignment_from_parts(
		{
			make_tuple( 2, 3, aln_1wzaA02_1zjaA02 ),
			make_tuple( 0, 3, aln_1g5aA03_1zjaA02 ),
			make_tuple( 0, 1, aln_1g5aA03_1r7aA02 )
		},
		all_proteins
	);

	// Check that the alignment matches the expected
	stringstream got_ss;
	write_alignment_as_fasta_alignment(
		got_ss,
		new_alignment,
		all_proteins
	);
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( got_ss, "got_ss", expected_aln_file );
}

BOOST_AUTO_TEST_SUITE_END()
