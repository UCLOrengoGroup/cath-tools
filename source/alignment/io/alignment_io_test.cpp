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
#include <boost/filesystem.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "alignment/io/alignment_io.hpp"
#include "chopping/region/region.hpp"
#include "common/boost_addenda/test/boost_check_equal_ranges.hpp"
#include "common/file/open_fstream.hpp"
#include "common/pair_insertion_operator.hpp"
#include "common/test_predicate/istreams_equal.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "ssap/ssap.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_io.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "test/global_test_constants.hpp"

#include <fstream>
#include <iostream>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace std;

//namespace std {
//	/// \brief Naughty addition of an insertion operator into std:: to get Boost.Test to output str_str_pairs
//	ostream & operator<<(ostream            &arg_os,          ///< ostream to which to output the str_str_pair
//	                     const str_str_pair &arg_str_str_pair ///< str_str_pair to output
//	                     ) {
//		arg_os << "pair[" << arg_str_str_pair.first  << ", ";
//		arg_os            << arg_str_str_pair.second << "]";
//		return arg_os;
//	}
//}

namespace cath {
	namespace test {

		/// \brief The alignment_io_test_suite_fixture to assist in testing alignment_io
		struct alignment_io_test_suite_fixture : protected global_test_constants {
		protected:
			~alignment_io_test_suite_fixture() noexcept = default;

		public:

			/// \brief TODOCUMENT
			void check_fasta_throws(const string &arg_fasta ///< The FASTA alignment string to parse
				                    ) {
				istringstream aln_ss( arg_fasta );
				BOOST_CHECK_THROW( read_ids_and_sequences_from_fasta( aln_ss ), runtime_error_exception );
			}

			/// \brief TODOCUMENT
			void check_fasta_gives_ids_and_seqs(const string           &arg_fasta,                 ///< The FASTA alignment string to parse
			                                    const str_str_pair_vec &arg_expcected_ids_and_seqs ///< The expected ids and sequences
			                                    ) {
				istringstream aln_ss( arg_fasta );
				const str_str_pair_vec got_ids_and_seqs = read_ids_and_sequences_from_fasta( aln_ss );
				BOOST_CHECK_EQUAL_RANGES( got_ids_and_seqs, arg_expcected_ids_and_seqs );
			}
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(alignment_io_test_suite, cath::test::alignment_io_test_suite_fixture)

/// \brief A sub test-suite for testing FASTA parsing
BOOST_AUTO_TEST_SUITE(fasta_test_suite)

/// \brief Check that a basic parse of FASTA input works as expected
BOOST_AUTO_TEST_CASE(fasta_parse_works) {
	check_fasta_gives_ids_and_seqs(
		">1d66B02\nTRAHLTEVESRLERL\n>1mkmA02\nGYKLIEYGSFVLRR-",
		{
			{ "1d66B02", "TRAHLTEVESRLERL" },
			{ "1mkmA02", "GYKLIEYGSFVLRR-" }
		}
	);
}

/// \brief Check that a FASTA parse correctly removes spaces and joins over newlines
BOOST_AUTO_TEST_CASE(fasta_parse_removes_spaces_and_newlines) {
	check_fasta_gives_ids_and_seqs(
		">1d66B02\nTRA\nH LTE V\nE SRLERL\n>1mkmA02\nGY  KL \nI EY\n GSFV\nLRR-",
		{
			{ "1d66B02", "TRAHLTEVESRLERL" },
			{ "1mkmA02", "GYKLIEYGSFVLRR-" }
		}
	);
}

/// \brief Check that a FASTA parse correctly upper-cases sequence letters
BOOST_AUTO_TEST_CASE(fasta_parse_uppercases_letters) {
	check_fasta_gives_ids_and_seqs(
		">1d66B02\nTRAHLteVESrLERL\n>1mkmA02\nGYklIEYGsfvlRR-",
		{
			{ "1d66B02", "TRAHLTEVESRLERL" },
			{ "1mkmA02", "GYKLIEYGSFVLRR-" }
		}
	);
}

/// \brief Check that a FASTA parse correctly throws a runtime_error_exception if any of the input contains non-printing characters
BOOST_AUTO_TEST_CASE(throws_on_non_printing_chars) {
	check_fasta_throws( string() + ">1d66B02\nTRAHLTEV" + "\x01\x05" + "ESRLERL\n>1mkmA02\nGYKLIEYGSFVL" + "\x0a\x15"  + "RR-" );
	check_fasta_throws( string() + ">1d66" + "\x01\x05" + "B02\nTRAHLTEVESRLERL\n>1mkm" + "\x0a\x15"  + "A02\nGYKLIEYGSFVLRR-" );
}

/// \brief Check that a FASTA parse correctly throws a runtime_error_exception if the first line doesn't start with a '>' symbol
BOOST_AUTO_TEST_CASE(throws_if_first_line_doesn_not_start_with_ge_symbol) {
	check_fasta_throws( "1d66B02\nTRAHLTEVESRLERL\n>1mkmA02\nGYKLIEYGSFVLRR-" );
}

/// \brief Check that a FASTA parse correctly throws a runtime_error_exception if any header line has an empty ID
BOOST_AUTO_TEST_CASE(throws_if_header_line_has_empty_id) {
	check_fasta_throws( ">1d66B02\nTRAHLTEVESRLERL\n>\nGYKLIEYGSFVLRR-" );
}

/// \brief Check that a FASTA parse correctly throws a runtime_error_exception if sequence lines contain characters other than spaces, letters or '-'
BOOST_AUTO_TEST_CASE(throws_if_sequence_line_contains_non_dash_or_letter_chars) {
	check_fasta_throws( ">1d66B02\nTRAHLTEVESRLERL\n>1mkmA02\nGYKLI1YGSFVLRR-" );
	check_fasta_throws( ">1d66B02\nTRAHLTEVESRLERL\n>1mkmA02\nGYKLI#YGSFVLRR-" );
	check_fasta_throws( ">1d66B02\nTRAHLTEVESRLERL\n>1mkmA02\nGYKLI$YGSFVLRR-" );
	check_fasta_throws( ">1d66B02\nTRAHLTEVESRLERL\n>1mkmA02\nGYKLI+YGSFVLRR-" );
	check_fasta_throws( ">1d66B02\nTRAHLTEVESRLERL\n>1mkmA02\nGYKLI@YGSFVLRR-" );
	check_fasta_throws( ">1d66B02\nTRAHLTEVESRLERL\n>1mkmA02\nGYKLI.YGSFVLRR-" );
}

BOOST_AUTO_TEST_SUITE_END()

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(alignment_legacy_input_output) {
	ostringstream err_ss;

	const protein protein_a = read_protein_from_dssp_and_pdb( EXAMPLE_A_DSSP_FILENAME(), EXAMPLE_A_PDB_FILENAME(), dssp_skip_policy::SKIP__BREAK_ANGLES, EXAMPLE_A_PDB_STEMNAME(), ostream_ref{ err_ss } );
	const protein protein_b = read_protein_from_dssp_and_pdb( EXAMPLE_B_DSSP_FILENAME(), EXAMPLE_B_PDB_FILENAME(), dssp_skip_policy::SKIP__BREAK_ANGLES, EXAMPLE_B_PDB_STEMNAME(), ostream_ref{ err_ss } );

	// Read the alignment file into a stringstream which can be used both as expected output
	// and as the istream for the parsing of the alignment
	stringstream expected_ss;
	ifstream alignment_file_stream;
	open_ifstream(alignment_file_stream, ALIGNMENT_FILE());
	expected_ss << alignment_file_stream.rdbuf();
	alignment_file_stream.close();

	// Parse the alignment from the stringstream, capturing stderr
	ostringstream test_stderr;
	const alignment my_aln = read_alignment_from_cath_ssap_legacy_format( expected_ss, protein_a, protein_b, ostream_ref{ test_stderr } );

	// Output the alignment to a stringstream
	ostringstream got_ss;
	output_alignment_to_cath_ssap_legacy_format(
		got_ss,
		my_aln,
		protein_a,
		protein_b
	);

	// Check that the data in the read+written alignment matches the original
	BOOST_CHECK_EQUAL( expected_ss.str(), got_ss.str() );
	BOOST_CHECK_EQUAL( err_ss.str(),      ""s          );
	
}

BOOST_AUTO_TEST_SUITE_END()
