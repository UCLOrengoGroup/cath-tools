/// \file
/// \brief The aligned_pair_score test suite

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

#include "aligned_pair_score.hpp"

#include <boost/filesystem/path.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "alignment/alignment.hpp"
#include "alignment/common_atom_selection_policy/common_atom_select_ca_policy.hpp"
#include "alignment/common_atom_selection_policy/common_atom_select_cb_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_min_score_policy.hpp"
#include "alignment/io/alignment_io.hpp"
#include "alignment/residue_score/residue_scorer.hpp"
#include "common/boost_addenda/log/log_to_ostream_guard.hpp"
#include "common/size_t_literal.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "score/aligned_pair_score/drmsd_score.hpp"
#include "score/aligned_pair_score/gsas_score.hpp"
#include "score/aligned_pair_score/lddt_score.hpp"
#include "score/aligned_pair_score/length_score.hpp"
#include "score/aligned_pair_score/mi_score.hpp"
#include "score/aligned_pair_score/overlap_score.hpp"
#include "score/aligned_pair_score/rmsd_score.hpp"
#include "score/aligned_pair_score/sas_score.hpp"
#include "score/aligned_pair_score/sequence_similarity_score.hpp"
#include "score/aligned_pair_score/si_score.hpp"
#include "score/aligned_pair_score/ssap_score.hpp"
#include "score/aligned_pair_score/substitution_matrix/blosum62_substitution_matrix.hpp"
#include "score/aligned_pair_score/substitution_matrix/match_substitution_matrix.hpp"
#include "score/aligned_pair_score/tm_score.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_list.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_list_factory.hpp"
#include "score/length_getter/length_of_first_getter.hpp"
#include "score/length_getter/length_of_longer_getter.hpp"
#include "score/length_getter/length_of_second_getter.hpp"
#include "score/length_getter/length_of_shorter_getter.hpp"
#include "score/length_getter/mean_length_getter.hpp"
#include "score/length_getter/num_aligned_length_getter.hpp"
#include "ssap/ssap.hpp"
#include "ssap/ssap_scores.hpp"
#include "structure/entry_querier/residue_querier.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/protein_source_file_set/protein_from_pdb.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "test/global_test_constants.hpp"
#include "test/test_tools.hpp"

#include <iterator>
#include <iostream>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::common::test;
using namespace cath::score;
using namespace std;

using boost::filesystem::path;

namespace cath {
	namespace test {

		/// \brief Fixture for the aligned_pair_score tests
		struct aligned_pair_score_fixture : protected global_test_constants {
		protected:
			~aligned_pair_score_fixture() noexcept = default;

		public:
			void check_ssap_scores(const alignment &,
			                       const protein &,
			                       const protein &,
			                       const score_value &,
			                       const score_value &,
			                       const score_value &) const;

			/// \brief TODOCUMENT
			const string NAME_1C55A               { "1c55A" };

			/// \brief TODOCUMENT
			const string NAME_1C56A               { "1c56A" };

			/// \brief TODOCUMENT
			const string NAME_1WMTA               { "1wmtA" };

			/// \brief TODOCUMENT
			const string NAME_1WT7A               { "1wt7A" };

			/// \brief TODOCUMENT
			const string NAME_1HYKA               { "1hykA" };

			/// \brief TODOCUMENT
			const path   ALIGN_PAIR_SCORE_TEST_DIR{ TEST_SOURCE_DATA_DIR() / "aligned_pair_score" };

			/// \brief TODOCUMENT
			const path   FASTA_ALN_1C55A_1C55A    { ALIGN_PAIR_SCORE_TEST_DIR / ( NAME_1C55A + "_" + NAME_1C55A + ".aln.fa" ) };

			/// \brief TODOCUMENT
			const path   FASTA_ALN_1C55A_1C56A    { ALIGN_PAIR_SCORE_TEST_DIR / ( NAME_1C55A + "_" + NAME_1C56A + ".aln.fa" ) };

			/// \brief TODOCUMENT
			const path   FASTA_ALN_1C55A_1WMTA    { ALIGN_PAIR_SCORE_TEST_DIR / ( NAME_1C55A + "_" + NAME_1WMTA + ".aln.fa" ) };

			/// \brief TODOCUMENT
			const path   FASTA_ALN_1C55A_1WT7A    { ALIGN_PAIR_SCORE_TEST_DIR / ( NAME_1C55A + "_" + NAME_1WT7A + ".aln.fa" ) };

			/// \brief TODOCUMENT
			const path   FASTA_ALN_1C55A_1HYKA    { ALIGN_PAIR_SCORE_TEST_DIR / ( NAME_1C55A + "_" + NAME_1HYKA + ".aln.fa" ) };

			ostringstream              parse_ss;
			const log_to_ostream_guard parse_log_guard      { parse_ss };

			const protein        protein_1c55A        { read_protein_from_files( protein_from_pdb(), ALIGN_PAIR_SCORE_TEST_DIR, "1c55A" ) };
			const protein        protein_1c56A        { read_protein_from_files( protein_from_pdb(), ALIGN_PAIR_SCORE_TEST_DIR, "1c56A" ) };
			const protein        protein_1hykA        { read_protein_from_files( protein_from_pdb(), ALIGN_PAIR_SCORE_TEST_DIR, "1hykA" ) };
			const protein        protein_1wmtA        { read_protein_from_files( protein_from_pdb(), ALIGN_PAIR_SCORE_TEST_DIR, "1wmtA" ) };
			const protein        protein_1wt7A        { read_protein_from_files( protein_from_pdb(), ALIGN_PAIR_SCORE_TEST_DIR, "1wt7A" ) };
			const protein_list   prot_list_1c55A_1c55A{ make_protein_list( { protein_1c55A, protein_1c55A } ) };
			const protein_list   prot_list_1c55A_1c56A{ make_protein_list( { protein_1c55A, protein_1c56A } ) };
			const protein_list   prot_list_1c55A_1hykA{ make_protein_list( { protein_1c55A, protein_1hykA } ) };
			const protein_list   prot_list_1c55A_1wmtA{ make_protein_list( { protein_1c55A, protein_1wmtA } ) };
			const protein_list   prot_list_1c55A_1wt7A{ make_protein_list( { protein_1c55A, protein_1wt7A } ) };
			const alignment      aln_1c55A_1c55A      { read_and_rescore_fasta_alignment( FASTA_ALN_1C55A_1C55A, prot_list_1c55A_1c55A, residue_scorer(), parse_ss ) };
			const alignment      aln_1c55A_1c56A      { read_and_rescore_fasta_alignment( FASTA_ALN_1C55A_1C56A, prot_list_1c55A_1c56A, residue_scorer(), parse_ss ) };
			const alignment      aln_1c55A_1hykA      { read_and_rescore_fasta_alignment( FASTA_ALN_1C55A_1HYKA, prot_list_1c55A_1hykA, residue_scorer(), parse_ss ) };
			const alignment      aln_1c55A_1wmtA      { read_and_rescore_fasta_alignment( FASTA_ALN_1C55A_1WMTA, prot_list_1c55A_1wmtA, residue_scorer(), parse_ss ) };
			const alignment      aln_1c55A_1wt7A      { read_and_rescore_fasta_alignment( FASTA_ALN_1C55A_1WT7A, prot_list_1c55A_1wt7A, residue_scorer(), parse_ss ) };
		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
void cath::test::aligned_pair_score_fixture::check_ssap_scores(const alignment   &arg_alignment,                      ///< TODOCUMENT
                                                               const protein     &arg_protein_a,                      ///< TODOCUMENT
                                                               const protein     &arg_protein_b,                      ///< TODOCUMENT
                                                               const score_value &arg_accurate_score_over_longer,     ///< TODOCUMENT
                                                               const score_value &arg_accurate_score_over_shorter,    ///< TODOCUMENT
                                                               const score_value &arg_accurate_score_over_num_aligned ///< TODOCUMENT
                                                               ) const {
	const ssap_scores the_ssap_scores = calculate_log_score( arg_alignment, arg_protein_a, arg_protein_b, residue_querier() );
	constexpr size_t exclude = 5;
	BOOST_CHECK_CLOSE(
		ssap_score( length_of_longer_getter(),   ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG,  ssap_score_accuracy::LOW, exclude ).calculate( arg_alignment, arg_protein_a, arg_protein_b ),
		the_ssap_scores.get_ssap_score_over_larger(),
		ACCURACY_PERCENTAGE()
	);
	BOOST_CHECK_CLOSE(
		ssap_score( length_of_shorter_getter(),  ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG,  ssap_score_accuracy::LOW, exclude ).calculate( arg_alignment, arg_protein_a, arg_protein_b ),
		the_ssap_scores.get_ssap_score_over_smaller(),
		ACCURACY_PERCENTAGE()
	);
	BOOST_CHECK_CLOSE(
		ssap_score( num_aligned_length_getter(), ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG,  ssap_score_accuracy::LOW, exclude ).calculate( arg_alignment, arg_protein_a, arg_protein_b ),
		the_ssap_scores.get_ssap_score_over_compared(),
		ACCURACY_PERCENTAGE()
	);

	BOOST_CHECK_CLOSE(
		ssap_score( length_of_longer_getter(),   ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG, ssap_score_accuracy::HIGH, exclude ).calculate( arg_alignment, arg_protein_a, arg_protein_b ),
		arg_accurate_score_over_longer,
		ACCURACY_PERCENTAGE()
	);
	BOOST_CHECK_CLOSE(
		ssap_score( length_of_shorter_getter(),  ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG, ssap_score_accuracy::HIGH, exclude ).calculate( arg_alignment, arg_protein_a, arg_protein_b ),
		arg_accurate_score_over_shorter,
		ACCURACY_PERCENTAGE()
	);
	BOOST_CHECK_CLOSE(
		ssap_score( num_aligned_length_getter(), ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG, ssap_score_accuracy::HIGH, exclude ).calculate( arg_alignment, arg_protein_a, arg_protein_b ),
		arg_accurate_score_over_num_aligned,
		ACCURACY_PERCENTAGE()
	);
}

/// \brief Test suite for testing various aligned pair scores
///
/// Here are the structures that are used in this test...
///
/// 1c55A versus 1c56A
/// ------------------
///  * highly similar (SSAP score ~97),
///  * identical sequences of 40 residues,
///  * identical residue numbering (1-40)
///
/// Results from TM-score webserver:
///
///     Structure1: A801003     Length=   40
///     Structure2: B801003     Length=   40 (by which all scores are normalized)
///     Number of residues in common=   40
///     RMSD of  the common residues=    0.605
///
///     TM-score    = 0.9132  (d0= 1.83)
///     MaxSub-score= 0.9723  (d0= 3.50)
///     GDT-TS-score= 0.9812 %(d<1)=0.9250 %(d<2)=1.0000 %(d<4)=1.0000 %(d<8)=1.0000
///     GDT-HA-score= 0.9063 %(d<0.5)=0.7000 %(d<1)=0.9250 %(d<2)=1.0000 %(d<4)=1.0000
///
/// 1c55A versus 1wt7A
/// ------------------
///  * structurally similar (SSAP score ~87),
///  * simple 1-1 alignment, but with 1wt7A having one extra residue at the end

///
/// 1c55A versus 1wmtA
/// ------------------
///  * structurally similar (SSAP score ~87),
///  * alignment with no gaps, but each sequence overhanging at one end
///
/// 1c55A versus 1hykA
/// ------------------
///  * similar?????????????,
///  * alignment with gaps in both sequences

BOOST_FIXTURE_TEST_SUITE(aligned_pair_score_test_suite, cath::test::aligned_pair_score_fixture)

/// \brief Sub test suite for values where the correct answer is known from some independent source
BOOST_AUTO_TEST_SUITE(known_correct_answers_suite)

BOOST_AUTO_TEST_SUITE(sequence_similarities_on_known)

/// \brief Check the sequence_similarities for 1c55A/1c55A are correct
BOOST_AUTO_TEST_CASE(seq_sims_1c55A_1c55A) {
	BOOST_CHECK_CLOSE( sequence_similarity_score(                          ).calculate( aln_1c55A_1c55A, protein_1c55A, protein_1c55A ), 100.0,          ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( sequence_similarity_score( make_subs_matrix_match() ).calculate( aln_1c55A_1c55A, protein_1c55A, protein_1c55A ), 100.0,          ACCURACY_PERCENTAGE()  );
}

/// \brief Check the sequence_similarities for 1c55A/1wmtA are correct
BOOST_AUTO_TEST_CASE(seq_sims_1c55A_1wmtA) {
	BOOST_CHECK_CLOSE( sequence_similarity_score(                          ).calculate( aln_1c55A_1wmtA, protein_1c55A, protein_1wmtA ),  35.0,          ACCURACY_PERCENTAGE()  ); // 14 out of 40
	BOOST_CHECK_CLOSE( sequence_similarity_score( make_subs_matrix_match() ).calculate( aln_1c55A_1wmtA, protein_1c55A, protein_1wmtA ), -20.0,          ACCURACY_PERCENTAGE()  ); // -8 out of 40
}

/// \brief Check the sequence_similarities for 1c55A/1wt7A are correct
BOOST_AUTO_TEST_CASE(seq_sims_1c55A_1wt7A) {
	BOOST_CHECK_CLOSE( sequence_similarity_score(                          ).calculate( aln_1c55A_1wt7A, protein_1c55A, protein_1wt7A ), 56.0975609756,  ACCURACY_PERCENTAGE()  ); // 23 out of 41
	BOOST_CHECK_CLOSE( sequence_similarity_score( make_subs_matrix_match() ).calculate( aln_1c55A_1wt7A, protein_1c55A, protein_1wt7A ), 14.6341463415,  ACCURACY_PERCENTAGE()  ); //  6 out of 41
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(overlaps_on_known)

/// \brief Check the overlaps for 1c55A/1c55A is correct
BOOST_AUTO_TEST_CASE(overlaps_1c55A_1c55A) {
	BOOST_CHECK_CLOSE( naive_overlap ().calculate( aln_1c55A_1c55A, protein_1c55A, protein_1c55A ), 100.0,          ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( local_overlap ().calculate( aln_1c55A_1c55A, protein_1c55A, protein_1c55A ), 100.0,          ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( global_overlap().calculate( aln_1c55A_1c55A, protein_1c55A, protein_1c55A ), 100.0,          ACCURACY_PERCENTAGE()  );
}

/// \brief Check the overlaps for 1c55A/1c56A is correct
BOOST_AUTO_TEST_CASE(overlaps_1c55A_1c56A) {
	BOOST_CHECK_CLOSE( naive_overlap ().calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 100.0,          ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( local_overlap ().calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 100.0,          ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( global_overlap().calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 100.0,          ACCURACY_PERCENTAGE()  );
}

/// \brief Check the overlaps for 1c55A/1wt7A is correct
BOOST_AUTO_TEST_CASE(overlaps_1c55A_1wt7A) {
	BOOST_CHECK_CLOSE( naive_overlap ().calculate( aln_1c55A_1wt7A, protein_1c55A, protein_1wt7A ),  97.5609756098, ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( local_overlap ().calculate( aln_1c55A_1wt7A, protein_1c55A, protein_1wt7A ), 100.0000000000, ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( global_overlap().calculate( aln_1c55A_1wt7A, protein_1c55A, protein_1wt7A ),  97.5609756098, ACCURACY_PERCENTAGE()  );
}

/// \brief Check the overlaps for 1c55A/1wmtA is correct
BOOST_AUTO_TEST_CASE(overlaps_1c55A_1wmtA) {
	BOOST_CHECK_CLOSE( naive_overlap ().calculate( aln_1c55A_1wmtA, protein_1c55A, protein_1wmtA ), 100.0,          ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( local_overlap ().calculate( aln_1c55A_1wmtA, protein_1c55A, protein_1wmtA ),  90.0,          ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( global_overlap().calculate( aln_1c55A_1wmtA, protein_1c55A, protein_1wmtA ),  90.0,          ACCURACY_PERCENTAGE()  );
}

/// \brief Check the overlaps for 1c55A/1hykA is correct
BOOST_AUTO_TEST_CASE(overlaps_1c55A_1hykA) {
//	cerr << precision (50)
	BOOST_CHECK_CLOSE( naive_overlap ().calculate( aln_1c55A_1hykA, protein_1c55A, protein_1hykA ),  86.9565217391, ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( local_overlap ().calculate( aln_1c55A_1hykA, protein_1c55A, protein_1hykA ),  77.5,          ACCURACY_PERCENTAGE()  );
	BOOST_CHECK_CLOSE( global_overlap().calculate( aln_1c55A_1hykA, protein_1c55A, protein_1hykA ),  67.3913043478, ACCURACY_PERCENTAGE()  );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(rmsd_on_known)

/// \brief Check the RMSD for 1c55A/1c55A is correct (TM-score webserver gave RMSD of 0.0)
BOOST_AUTO_TEST_CASE(rmsd_1c55A_1c55A) {
	BOOST_CHECK_SMALL( rmsd_score().calculate( aln_1c55A_1c55A, protein_1c55A, protein_1c55A ),                     ACCURACY_PERCENTAGE()  );
}

/// \brief Check the RMSD for 1c55A/1c56A is correct (TM-score webserver gave RMSD of 0.605)
BOOST_AUTO_TEST_CASE(rmsd_1c55A_1c56A) {
	BOOST_CHECK_CLOSE( rmsd_score().calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.60524408098498938546327963194926269352436065673828,  ACCURACY_PERCENTAGE() );
}

/// \brief Check the RMSD for 1c55A/1wt7A is correct (TM-score webserver gave RMSD of 1.99)
BOOST_AUTO_TEST_CASE(rmsd_1c55A_1wt7A) {
	BOOST_CHECK_CLOSE( rmsd_score().calculate( aln_1c55A_1wt7A, protein_1c55A, protein_1wt7A ), 4.9863504664588501, ACCURACY_PERCENTAGE() );
}

/// \brief Check the RMSD for 1c55A/1wmtA is correct (TM-score webserver gave RMSD of 1.76)
BOOST_AUTO_TEST_CASE(rmsd_1c55A_1wmtA) {
	BOOST_CHECK_CLOSE( rmsd_score().calculate( aln_1c55A_1wmtA, protein_1c55A, protein_1wmtA ), 2.040888817288554, ACCURACY_PERCENTAGE() );
}

/// \brief Check the RMSD for 1c55A/1hykA is correct (TM-score webserver gave RMSD of 0.605)
BOOST_AUTO_TEST_CASE(rmsd_1c55A_1hykA) {
	BOOST_CHECK_CLOSE( rmsd_score().calculate( aln_1c55A_1hykA, protein_1c55A, protein_1hykA ), 5.1230997129218308, ACCURACY_PERCENTAGE() );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(tm_score_on_known)

/// \brief Check the TM-score for 1c55A/1c55A is correct
BOOST_AUTO_TEST_CASE(tm_score_1c55A_1c55A) {
	BOOST_CHECK_EQUAL( tm_score().calculate  ( aln_1c55A_1c55A, protein_1c55A, protein_1c55A ), 1.0 );
}

/// \brief Check the TM-score for 1c55A/1c56A is correct
BOOST_AUTO_TEST_CASE(tm_score_1c55A_1c56A) {
	BOOST_CHECK_CLOSE( tm_score().calculate  ( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.9131980677380123, ACCURACY_PERCENTAGE() );
}

/// \brief Check the TM-score for 1c55A/1wt7A is correct
BOOST_AUTO_TEST_CASE(tm_score_1c55A_1wt7A) {
	BOOST_WARN_CLOSE( tm_score().calculate  ( aln_1c55A_1wt7A, protein_1c55A, protein_1wt7A ), 0.58304,            ACCURACY_PERCENTAGE() );
	BOOST_CHECK( true );
}

/// \brief Check the TM-score for 1c55A/1wmtA is correct
BOOST_AUTO_TEST_CASE(tm_score_1c55A_1wmtA) {
	BOOST_WARN_CLOSE( tm_score().calculate  ( aln_1c55A_1wmtA, protein_1c55A, protein_1wmtA ), 0.54228,            ACCURACY_PERCENTAGE() );
	BOOST_CHECK( true );
}

/// \brief Check the TM-score for 1c55A/1hykA is correct
BOOST_AUTO_TEST_CASE(tm_score_1c55A_1hykA) {
	BOOST_WARN_CLOSE( tm_score().calculate  ( aln_1c55A_1hykA, protein_1c55A, protein_1hykA ), 0.34689,            ACCURACY_PERCENTAGE() );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()



/// \brief Sub test suite where the tests are just regression tests for which the correct answer has not been taken from any independent source
BOOST_AUTO_TEST_SUITE(regression_suite)

/// \brief Check that the
BOOST_AUTO_TEST_CASE(drmsd_score_1c55A_1c56A) {
	BOOST_CHECK_CLOSE( drmsd_score(                                                                                   ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.452877184244952935011, ACCURACY_PERCENTAGE() );
	BOOST_CHECK_CLOSE( drmsd_score( common_residue_select_all_policy(),                common_atom_select_ca_policy() ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.452877184244952935011, ACCURACY_PERCENTAGE() );
	BOOST_CHECK_CLOSE( drmsd_score( common_residue_select_all_policy(),                common_atom_select_cb_policy() ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.506673194154965700342, ACCURACY_PERCENTAGE() );
	BOOST_CHECK_CLOSE( drmsd_score( common_residue_select_best_score_percent_policy(), common_atom_select_ca_policy() ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.258779215015056851534, ACCURACY_PERCENTAGE() );
	BOOST_CHECK_CLOSE( drmsd_score( common_residue_select_best_score_percent_policy(), common_atom_select_cb_policy() ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.291026502291049260496, ACCURACY_PERCENTAGE() );
	BOOST_CHECK_CLOSE( drmsd_score( common_residue_select_min_score_policy(),          common_atom_select_ca_policy() ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.452877184244952935011, ACCURACY_PERCENTAGE() );
	BOOST_CHECK_CLOSE( drmsd_score( common_residue_select_min_score_policy(),          common_atom_select_cb_policy() ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.506673194154965700342, ACCURACY_PERCENTAGE() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(gsas_score_1c55A_1c56A) {
	BOOST_CHECK_CLOSE(  gsas_score().calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 1.513110202462473408147, ACCURACY_PERCENTAGE() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(lddt_score_1c55A_1c56A) {
	BOOST_CHECK_CLOSE(  lddt_score().calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.944196428571428603149, ACCURACY_PERCENTAGE() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(mi_score_1c55A_1c56A) {
	BOOST_CHECK_CLOSE(    mi_score().calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.287493543599852463544, ACCURACY_PERCENTAGE() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(sas_score_1c55A_1c56A) {
	BOOST_CHECK_CLOSE(   sas_score().calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 1.513110202462473408147, ACCURACY_PERCENTAGE() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(si_score_1c55A_1c56A) {
	BOOST_CHECK_CLOSE(    si_score().calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 0.605244080984989385463, ACCURACY_PERCENTAGE() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(length_score_1c55A_1c56A) {
	BOOST_CHECK_EQUAL( length_score( length_of_first_getter()    ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 40 );
	BOOST_CHECK_EQUAL( length_score( length_of_longer_getter()   ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 40 );
	BOOST_CHECK_EQUAL( length_score( length_of_second_getter()   ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 40 );
	BOOST_CHECK_EQUAL( length_score( length_of_shorter_getter()  ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 40 );
	BOOST_CHECK_EQUAL( length_score( num_aligned_length_getter() ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 40 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(sequence_identity_1c55A_1c56A) {
	BOOST_CHECK_EQUAL( sequence_similarity_score(                             ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 100.0000000000 );
	BOOST_CHECK_EQUAL( sequence_similarity_score( make_subs_matrix_blosum62() ).calculate( aln_1c55A_1c56A, protein_1c55A, protein_1c56A ), 620.0 / 11.0 );

	const ssap_scores the_ssap_scores = calculate_log_score( aln_1c55A_1c56A, protein_1c55A, protein_1c56A, residue_querier() );
	BOOST_CHECK_EQUAL( the_ssap_scores.get_num_aligned_pairs(),                     40_z );
	BOOST_CHECK_EQUAL( the_ssap_scores.get_seq_id(),                               100.0 );
	BOOST_CHECK_EQUAL( the_ssap_scores.get_percentage_aligned_pairs_over_larger(), 100.0 );
}

/// \brief Test SSAP scores for 1c55A versus 1c56A (very similar, identical sequences, identical numbering, identity alignment)
BOOST_AUTO_TEST_CASE(ssap_score_score_1c55A_1c56A) {
	check_ssap_scores(
		aln_1c55A_1c56A,
		protein_1c55A,
		protein_1c56A,
		97.563912614059461,
		97.563912614059461,
		97.563912614059461
	);
}

/// \brief Test SSAP scores for 1c55A versus 1wt7A (similar, non-gapped alignment but one sequence has an extra residue at the end)
BOOST_AUTO_TEST_CASE(ssap_score_score_1c55A_1wt7A) {
	check_ssap_scores(
		aln_1c55A_1wt7A,
		protein_1c55A,
		protein_1wt7A,
		87.467595428797367,
		87.99587255553206,
		87.99587255553206
	);
}

/// \brief Test SSAP scores for 1c55A versus 1wmtA (similar, non-gapped alignment but each sequence overhangs by a few residues at one end)
BOOST_AUTO_TEST_CASE(ssap_score_score_1c55A_1wmtA) {
	check_ssap_scores(
		aln_1c55A_1wmtA,
		protein_1c55A,
		protein_1wmtA,
		87.049696689395191,
		87.049696689395191,
		89.328153728699561
	);
}

/// \brief Test SSAP scores for 1c55A versus 1hykA (similar, alignment has gaps in both structures)
BOOST_AUTO_TEST_CASE(ssap_score_score_1c55A_1hykA) {
	check_ssap_scores(
		aln_1c55A_1hykA,
		protein_1c55A,
		protein_1hykA,
		71.474643160901948,
		74.439056352018923,
		80.028230245602572
	);
}

/// \brief Test that several key human_friendly_short_names haven't changed from these expected values
BOOST_AUTO_TEST_CASE(human_friendly_short_name_regression) {
	BOOST_CHECK_EQUAL(  drmsd_score(                                       common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "dRMSD"                                                        );
	BOOST_CHECK_EQUAL(  drmsd_score(                                       common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "dRMSD.cb_atoms"                                               );
	BOOST_CHECK_EQUAL(  drmsd_score(                                       common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "dRMSD.select_min_score[-0.25]"                                );
	BOOST_CHECK_EQUAL(  drmsd_score(                                       common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "dRMSD.select_min_score[-0.25].cb_atoms"                       );

	BOOST_CHECK_EQUAL(   gsas_score(                                       common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "GSAS"                                                         );
	BOOST_CHECK_EQUAL(   gsas_score(                                       common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "GSAS.cb_atoms"                                                );
	BOOST_CHECK_EQUAL(   gsas_score(                                       common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "GSAS.select_min_score[-0.25]"                                 );
	BOOST_CHECK_EQUAL(   gsas_score(                                       common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "GSAS.select_min_score[-0.25].cb_atoms"                        );

	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::DEFAULT_AVG, common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "lDDT.threshold[STD_MEAN]"                                     );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::DEFAULT_AVG, common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "lDDT.threshold[STD_MEAN].cb_atoms"                            );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::DEFAULT_AVG, common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "lDDT.threshold[STD_MEAN].select_min_score[-0.25]"             );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::DEFAULT_AVG, common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "lDDT.threshold[STD_MEAN].select_min_score[-0.25].cb_atoms"    );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::HALF_A,      common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "lDDT.threshold[0.5]"                                          );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::HALF_A,      common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "lDDT.threshold[0.5].cb_atoms"                                 );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::HALF_A,      common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "lDDT.threshold[0.5].select_min_score[-0.25]"                  );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::HALF_A,      common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "lDDT.threshold[0.5].select_min_score[-0.25].cb_atoms"         );

	BOOST_CHECK_EQUAL( length_score( length_of_longer_getter()                                                                                      ).human_friendly_short_name(), "longer_protein_length"                                        );
	BOOST_CHECK_EQUAL( length_score( length_of_shorter_getter()                                                                                     ).human_friendly_short_name(), "shorter_protein_length"                                       );
	BOOST_CHECK_EQUAL( length_score( mean_length_getter()                                                                                           ).human_friendly_short_name(), "mean_protein_length"                                          );

	BOOST_CHECK_EQUAL(     mi_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "MIMAX"                                                        );
	BOOST_CHECK_EQUAL(     mi_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "MIMAX.cb_atoms"                                               );
	BOOST_CHECK_EQUAL(     mi_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "MIMAX.select_min_score[-0.25]"                                );
	BOOST_CHECK_EQUAL(     mi_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "MIMAX.select_min_score[-0.25].cb_atoms"                       );
	BOOST_CHECK_EQUAL(     mi_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "MI"                                    );
	BOOST_CHECK_EQUAL(     mi_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "MI.cb_atoms"                           );
	BOOST_CHECK_EQUAL(     mi_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "MI.select_min_score[-0.25]"            );
	BOOST_CHECK_EQUAL(     mi_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "MI.select_min_score[-0.25].cb_atoms"   );
	BOOST_CHECK_EQUAL(     mi_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "MI.mean_protein_length"                                       );
	BOOST_CHECK_EQUAL(     mi_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "MI.mean_protein_length.cb_atoms"                              );
	BOOST_CHECK_EQUAL(     mi_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "MI.mean_protein_length.select_min_score[-0.25]"               );
	BOOST_CHECK_EQUAL(     mi_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "MI.mean_protein_length.select_min_score[-0.25].cb_atoms"      );

	BOOST_CHECK_EQUAL(   rmsd_score(                                       common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "RMSD"                                                         );
	BOOST_CHECK_EQUAL(   rmsd_score(                                       common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "RMSD.cb_atoms"                                                );
	BOOST_CHECK_EQUAL(   rmsd_score(                                       common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "RMSD.select_min_score[-0.25]"                                 );
	BOOST_CHECK_EQUAL(   rmsd_score(                                       common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "RMSD.select_min_score[-0.25].cb_atoms"                        );

	BOOST_CHECK_EQUAL(     si_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "SIMAX"                                                        );
	BOOST_CHECK_EQUAL(     si_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "SIMAX.cb_atoms"                                               );
	BOOST_CHECK_EQUAL(     si_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "SIMAX.select_min_score[-0.25]"                                );
	BOOST_CHECK_EQUAL(     si_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "SIMAX.select_min_score[-0.25].cb_atoms"                       );
	BOOST_CHECK_EQUAL(     si_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "SI"                                    );
	BOOST_CHECK_EQUAL(     si_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "SI.cb_atoms"                           );
	BOOST_CHECK_EQUAL(     si_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "SI.select_min_score[-0.25]"            );
	BOOST_CHECK_EQUAL(     si_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "SI.select_min_score[-0.25].cb_atoms"   );
	BOOST_CHECK_EQUAL(     si_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "SI.mean_protein_length"                                       );
	BOOST_CHECK_EQUAL(     si_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "SI.mean_protein_length.cb_atoms"                              );
	BOOST_CHECK_EQUAL(     si_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "SI.mean_protein_length.select_min_score[-0.25]"               );
	BOOST_CHECK_EQUAL(     si_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "SI.mean_protein_length.select_min_score[-0.25].cb_atoms"      );

	BOOST_CHECK_EQUAL(    sas_score(                                       common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "SAS"                                                          );
	BOOST_CHECK_EQUAL(    sas_score(                                       common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "SAS.cb_atoms"                                                 );
	BOOST_CHECK_EQUAL(    sas_score(                                       common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "SAS.select_min_score[-0.25]"                                  );
	BOOST_CHECK_EQUAL(    sas_score(                                       common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "SAS.select_min_score[-0.25].cb_atoms"                         );

	BOOST_CHECK_EQUAL(   ssap_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "ssap"                                                         );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "ssap.cb_atoms"                                                );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "ssap.select_min_score[-0.25]"                                 );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "ssap.select_min_score[-0.25].cb_atoms"                        );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "ssap.shorter_protein_length"                                  );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "ssap.shorter_protein_length.cb_atoms"                         );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "ssap.shorter_protein_length.select_min_score[-0.25]"          );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "ssap.shorter_protein_length.select_min_score[-0.25].cb_atoms" );
	BOOST_CHECK_EQUAL(   ssap_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_ca_policy() ).human_friendly_short_name(), "ssap.mean_protein_length"                                     );
	BOOST_CHECK_EQUAL(   ssap_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_cb_policy() ).human_friendly_short_name(), "ssap.mean_protein_length.cb_atoms"                            );
	BOOST_CHECK_EQUAL(   ssap_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).human_friendly_short_name(), "ssap.mean_protein_length.select_min_score[-0.25]"             );
	BOOST_CHECK_EQUAL(   ssap_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).human_friendly_short_name(), "ssap.mean_protein_length.select_min_score[-0.25].cb_atoms"    );

	BOOST_CHECK_EQUAL(  sequence_similarity_score( length_of_longer_getter(),      make_subs_matrix_blosum62()                                      ).human_friendly_short_name(), "sequence_id.blosum62"                                         );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( length_of_longer_getter(),      make_subs_matrix_identity()                                      ).human_friendly_short_name(), "sequence_id.identity"                                         );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( length_of_shorter_getter(),     make_subs_matrix_blosum62()                                      ).human_friendly_short_name(), "sequence_id.blosum62.shorter_protein_length"                  );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( length_of_shorter_getter(),     make_subs_matrix_identity()                                      ).human_friendly_short_name(), "sequence_id.identity.shorter_protein_length"                  );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( mean_length_getter(),           make_subs_matrix_blosum62()                                      ).human_friendly_short_name(), "sequence_id.blosum62.mean_protein_length"                     );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( mean_length_getter(),           make_subs_matrix_identity()                                      ).human_friendly_short_name(), "sequence_id.identity.mean_protein_length"                     );

	BOOST_CHECK_EQUAL( naive_overlap ( ).human_friendly_short_name(), "overlap.shorter_protein_length.longer_protein_length" );
	BOOST_CHECK_EQUAL( local_overlap ( ).human_friendly_short_name(), "overlap.num_aligned_residues.shorter_protein_length"  );
	BOOST_CHECK_EQUAL( global_overlap( ).human_friendly_short_name(), "overlap.num_aligned_residues.longer_protein_length"   );
}


/// \brief Test that several key short_names haven't changed from these expected values
BOOST_AUTO_TEST_CASE(full_short_name_regression) {
	BOOST_CHECK_EQUAL(  drmsd_score(                                       common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "dRMSD.select_all.ca_atoms"                                    );
	BOOST_CHECK_EQUAL(  drmsd_score(                                       common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "dRMSD.select_all.cb_atoms"                                    );
	BOOST_CHECK_EQUAL(  drmsd_score(                                       common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "dRMSD.select_min_score[-0.25].ca_atoms"                       );
	BOOST_CHECK_EQUAL(  drmsd_score(                                       common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "dRMSD.select_min_score[-0.25].cb_atoms"                       );

	BOOST_CHECK_EQUAL(   gsas_score(                                       common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "GSAS.select_all.ca_atoms"                                     );
	BOOST_CHECK_EQUAL(   gsas_score(                                       common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "GSAS.select_all.cb_atoms"                                     );
	BOOST_CHECK_EQUAL(   gsas_score(                                       common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "GSAS.select_min_score[-0.25].ca_atoms"                        );
	BOOST_CHECK_EQUAL(   gsas_score(                                       common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "GSAS.select_min_score[-0.25].cb_atoms"                        );

	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::DEFAULT_AVG, common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "lDDT.threshold[STD_MEAN].select_all.ca_atoms"                 );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::DEFAULT_AVG, common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "lDDT.threshold[STD_MEAN].select_all.cb_atoms"                 );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::DEFAULT_AVG, common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "lDDT.threshold[STD_MEAN].select_min_score[-0.25].ca_atoms"    );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::DEFAULT_AVG, common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "lDDT.threshold[STD_MEAN].select_min_score[-0.25].cb_atoms"    );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::HALF_A,      common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "lDDT.threshold[0.5].select_all.ca_atoms"                      );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::HALF_A,      common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "lDDT.threshold[0.5].select_all.cb_atoms"                      );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::HALF_A,      common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "lDDT.threshold[0.5].select_min_score[-0.25].ca_atoms"         );
	BOOST_CHECK_EQUAL(   lddt_score( lddt_distance_threshold::HALF_A,      common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "lDDT.threshold[0.5].select_min_score[-0.25].cb_atoms"         );

	BOOST_CHECK_EQUAL( length_score( length_of_longer_getter()                                                                                      ).full_short_name(), "longer_protein_length"                                        );
	BOOST_CHECK_EQUAL( length_score( length_of_shorter_getter()                                                                                     ).full_short_name(), "shorter_protein_length"                                       );
	BOOST_CHECK_EQUAL( length_score( mean_length_getter()                                                                                           ).full_short_name(), "mean_protein_length"                                          );

	BOOST_CHECK_EQUAL(     mi_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "MIMAX.longer_protein_length.select_all.ca_atoms"              );
	BOOST_CHECK_EQUAL(     mi_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "MIMAX.longer_protein_length.select_all.cb_atoms"              );
	BOOST_CHECK_EQUAL(     mi_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "MIMAX.longer_protein_length.select_min_score[-0.25].ca_atoms" );
	BOOST_CHECK_EQUAL(     mi_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "MIMAX.longer_protein_length.select_min_score[-0.25].cb_atoms" );
	BOOST_CHECK_EQUAL(     mi_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "MI.shorter_protein_length.select_all.ca_atoms"                );
	BOOST_CHECK_EQUAL(     mi_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "MI.shorter_protein_length.select_all.cb_atoms"                );
	BOOST_CHECK_EQUAL(     mi_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "MI.shorter_protein_length.select_min_score[-0.25].ca_atoms"   );
	BOOST_CHECK_EQUAL(     mi_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "MI.shorter_protein_length.select_min_score[-0.25].cb_atoms"   );
	BOOST_CHECK_EQUAL(     mi_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "MI.mean_protein_length.select_all.ca_atoms"                   );
	BOOST_CHECK_EQUAL(     mi_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "MI.mean_protein_length.select_all.cb_atoms"                   );
	BOOST_CHECK_EQUAL(     mi_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "MI.mean_protein_length.select_min_score[-0.25].ca_atoms"      );
	BOOST_CHECK_EQUAL(     mi_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "MI.mean_protein_length.select_min_score[-0.25].cb_atoms"      );

	BOOST_CHECK_EQUAL(   rmsd_score(                                       common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "RMSD.select_all.ca_atoms"                                     );
	BOOST_CHECK_EQUAL(   rmsd_score(                                       common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "RMSD.select_all.cb_atoms"                                     );
	BOOST_CHECK_EQUAL(   rmsd_score(                                       common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "RMSD.select_min_score[-0.25].ca_atoms"                        );
	BOOST_CHECK_EQUAL(   rmsd_score(                                       common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "RMSD.select_min_score[-0.25].cb_atoms"                        );

	BOOST_CHECK_EQUAL(     si_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "SIMAX.longer_protein_length.select_all.ca_atoms"              );
	BOOST_CHECK_EQUAL(     si_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "SIMAX.longer_protein_length.select_all.cb_atoms"              );
	BOOST_CHECK_EQUAL(     si_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "SIMAX.longer_protein_length.select_min_score[-0.25].ca_atoms" );
	BOOST_CHECK_EQUAL(     si_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "SIMAX.longer_protein_length.select_min_score[-0.25].cb_atoms" );
	BOOST_CHECK_EQUAL(     si_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "SI.shorter_protein_length.select_all.ca_atoms"                );
	BOOST_CHECK_EQUAL(     si_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "SI.shorter_protein_length.select_all.cb_atoms"                );
	BOOST_CHECK_EQUAL(     si_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "SI.shorter_protein_length.select_min_score[-0.25].ca_atoms"   );
	BOOST_CHECK_EQUAL(     si_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "SI.shorter_protein_length.select_min_score[-0.25].cb_atoms"   );
	BOOST_CHECK_EQUAL(     si_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "SI.mean_protein_length.select_all.ca_atoms"                   );
	BOOST_CHECK_EQUAL(     si_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "SI.mean_protein_length.select_all.cb_atoms"                   );
	BOOST_CHECK_EQUAL(     si_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "SI.mean_protein_length.select_min_score[-0.25].ca_atoms"      );
	BOOST_CHECK_EQUAL(     si_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "SI.mean_protein_length.select_min_score[-0.25].cb_atoms"      );

	BOOST_CHECK_EQUAL(    sas_score(                                       common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "SAS.select_all.ca_atoms"                                      );
	BOOST_CHECK_EQUAL(    sas_score(                                       common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "SAS.select_all.cb_atoms"                                      );
	BOOST_CHECK_EQUAL(    sas_score(                                       common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "SAS.select_min_score[-0.25].ca_atoms"                         );
	BOOST_CHECK_EQUAL(    sas_score(                                       common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "SAS.select_min_score[-0.25].cb_atoms"                         );

	BOOST_CHECK_EQUAL(   ssap_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "ssap.longer_protein_length.select_all.ca_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"               );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_longer_getter(),            common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "ssap.longer_protein_length.select_all.cb_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"               );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "ssap.longer_protein_length.select_min_score[-0.25].ca_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"  );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_longer_getter(),            common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "ssap.longer_protein_length.select_min_score[-0.25].cb_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"  );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "ssap.shorter_protein_length.select_all.ca_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"              );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_shorter_getter(),           common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "ssap.shorter_protein_length.select_all.cb_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"              );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "ssap.shorter_protein_length.select_min_score[-0.25].ca_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code" );
	BOOST_CHECK_EQUAL(   ssap_score( length_of_shorter_getter(),           common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "ssap.shorter_protein_length.select_min_score[-0.25].cb_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code" );
	BOOST_CHECK_EQUAL(   ssap_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_ca_policy() ).full_short_name(), "ssap.mean_protein_length.select_all.ca_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"                 );
	BOOST_CHECK_EQUAL(   ssap_score( mean_length_getter(),                 common_residue_select_all_policy(),       common_atom_select_cb_policy() ).full_short_name(), "ssap.mean_protein_length.select_all.cb_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"                 );
	BOOST_CHECK_EQUAL(   ssap_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_ca_policy() ).full_short_name(), "ssap.mean_protein_length.select_min_score[-0.25].ca_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"    );
	BOOST_CHECK_EQUAL(   ssap_score( mean_length_getter(),                 common_residue_select_min_score_policy(), common_atom_select_cb_policy() ).full_short_name(), "ssap.mean_protein_length.select_min_score[-0.25].cb_atoms.ssap_score_post_processing::complex_normalise_then_log.low_accuracy.num_excluded_on_sides:5.distance_score_formula_used_in_previous_code"    );

	BOOST_CHECK_EQUAL( naive_overlap ( ).full_short_name(), "overlap.shorter_protein_length.longer_protein_length"                    );
	BOOST_CHECK_EQUAL( local_overlap ( ).full_short_name(), "overlap.num_aligned_residues.select_all.ca_atoms.shorter_protein_length" );
	BOOST_CHECK_EQUAL( global_overlap( ).full_short_name(), "overlap.num_aligned_residues.select_all.ca_atoms.longer_protein_length"  );

	BOOST_CHECK_EQUAL(  sequence_similarity_score( length_of_longer_getter(),      make_subs_matrix_blosum62()                                      ).full_short_name(), "sequence_id.blosum62.longer_protein_length"                   );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( length_of_longer_getter(),      make_subs_matrix_identity()                                      ).full_short_name(), "sequence_id.identity.longer_protein_length"                   );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( length_of_shorter_getter(),     make_subs_matrix_blosum62()                                      ).full_short_name(), "sequence_id.blosum62.shorter_protein_length"                  );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( length_of_shorter_getter(),     make_subs_matrix_identity()                                      ).full_short_name(), "sequence_id.identity.shorter_protein_length"                  );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( mean_length_getter(),           make_subs_matrix_blosum62()                                      ).full_short_name(), "sequence_id.blosum62.mean_protein_length"                     );
	BOOST_CHECK_EQUAL(  sequence_similarity_score( mean_length_getter(),           make_subs_matrix_identity()                                      ).full_short_name(), "sequence_id.identity.mean_protein_length"                     );
}

BOOST_AUTO_TEST_SUITE_END()

/// \brief Sub test suite for checks of functionality
BOOST_AUTO_TEST_SUITE(functionality)

/// \brief Test equality/inequality for rmsd_score pair that previously failed because it didn't compare its score_common_coord_handler
BOOST_AUTO_TEST_CASE(check_equality_for_rmsd) {
	check_equality_operators_on_diff_vals(
		rmsd_score(),
		rmsd_score( common_residue_select_all_policy(), common_atom_select_cb_policy() )
	);
}

/// \brief Test equality/inequality for lddt_score pair that previously failed because it didn't compare its score_common_coord_handler
BOOST_AUTO_TEST_CASE(check_equality_for_lddt_score) {
	check_equality_operators_on_diff_vals(
		lddt_score( lddt_distance_threshold::DEFAULT_AVG ),
		lddt_score( lddt_distance_threshold::DEFAULT_AVG, common_residue_select_all_policy(), common_atom_select_cb_policy() )
	);
}

/// \brief Test that equality/inequality for pairs of members in standard aligned_pair_score_lists
BOOST_AUTO_TEST_CASE(check_equality_over_lists) {
	// Loop over three lists of scores, each of which should only contain unique aligned_pair_scores
	for (const aligned_pair_score_list &the_scores : {
//	                                                   make_full_aligned_pair_score_list(), // This is commented out because it's a bit too slow
//	                                                   make_old_ssap_aligned_pair_score_list(),
	                                                   make_default_aligned_pair_score_list() } ) {
		// Grab the number of scores in this list and require that it's a decent number
		const size_t num_scores = the_scores.size();
		BOOST_REQUIRE_GE( num_scores, 6 );

		// Loop over the possible non-equal pairs of score indices
		check_equality_operators_on_diff_vals_range( the_scores );
	}
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

