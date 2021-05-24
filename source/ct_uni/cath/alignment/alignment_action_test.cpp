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

#include "alignment_action.hpp"

#include <filesystem>

#include <boost/test/unit_test.hpp>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/io/alignment_io.hpp"
#include "cath/file/pdb/dssp_skip_policy.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_io.hpp"
#include "cath/structure/protein/protein_list.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"
#include "cath/test/global_test_constants.hpp"
#include "cath/test/predicate/istream_and_file_equal.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::file;
using namespace ::std;

using ::std::filesystem::path;

namespace {

	/// \brief The alignment_action_test_suite_fixture to assist in testing alignment_action
	struct alignment_action_test_suite_fixture : protected global_test_constants {
	private:
		ostringstream test_stderr;

	protected:
		~alignment_action_test_suite_fixture() noexcept;

		const path         root_dir           = { TEST_MULTI_SSAP_SUPERPOSE_DIR() };
		const path         expected_aln_file  = { root_dir / "expected_glue_output.fa" };
		const protein      protein_1g5aA03    = { read_protein_from_dssp_and_pdb( root_dir / "1g5aA03.dssp", root_dir / "1g5aA03", dssp_skip_policy::SKIP__BREAK_ANGLES, "1g5aA03", ostream_ref( test_stderr ) ) };
		const protein      protein_1r7aA02    = { read_protein_from_dssp_and_pdb( root_dir / "1r7aA02.dssp", root_dir / "1r7aA02", dssp_skip_policy::SKIP__BREAK_ANGLES, "1r7aA02", ostream_ref( test_stderr ) ) };
		const protein      protein_1wzaA02    = { read_protein_from_dssp_and_pdb( root_dir / "1wzaA02.dssp", root_dir / "1wzaA02", dssp_skip_policy::SKIP__BREAK_ANGLES, "1wzaA02", ostream_ref( test_stderr ) ) };
		const protein      protein_1zjaA02    = { read_protein_from_dssp_and_pdb( root_dir / "1zjaA02.dssp", root_dir / "1zjaA02", dssp_skip_policy::SKIP__BREAK_ANGLES, "1zjaA02", ostream_ref( test_stderr ) ) };
		const protein_list all_proteins       = { make_protein_list( { protein_1g5aA03, protein_1r7aA02, protein_1wzaA02, protein_1zjaA02 } ) };
		const alignment    aln_1wzaA02_1zjaA02= { read_alignment_from_cath_ssap_legacy_format( root_dir / "1wzaA021zjaA02.list", protein_1wzaA02, protein_1zjaA02, ostream_ref( test_stderr ) ) };
		const alignment    aln_1g5aA03_1zjaA02= { read_alignment_from_cath_ssap_legacy_format( root_dir / "1g5aA031zjaA02.list", protein_1g5aA03, protein_1zjaA02, ostream_ref( test_stderr ) ) };
		const alignment    aln_1g5aA03_1r7aA02= { read_alignment_from_cath_ssap_legacy_format( root_dir / "1g5aA031r7aA02.list", protein_1g5aA03, protein_1r7aA02, ostream_ref( test_stderr ) ) };
	};

	/// \brief TODOCUMENT
	alignment_action_test_suite_fixture::~alignment_action_test_suite_fixture() noexcept {
		try {
			BOOST_CHECK_EQUAL( test_stderr.str(), ""s );
		}
		catch (...) { /// Prevent the destructor throwing any exceptions
		}
	}

} // namespace

BOOST_FIXTURE_TEST_SUITE(alignment_action_test_suite, alignment_action_test_suite_fixture)

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
		all_proteins,
		aln_glue_style::SIMPLY
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
