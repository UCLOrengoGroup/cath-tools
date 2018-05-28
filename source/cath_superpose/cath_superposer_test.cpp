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

#include "cath_superposer.hpp"

#include <boost/log/trivial.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "cath_superpose/options/cath_superpose_options.hpp"
#include "chopping/domain/domain.hpp"
#include "common/argc_argv_faker.hpp"
#include "common/boost_addenda/log/stringstream_log_sink.hpp"
#include "common/file/open_fstream.hpp"
#include "common/file/temp_file.hpp"
#include "common/regex/regex_replace_file.hpp"
#include "test/global_test_constants.hpp"
#include "test/predicate/files_equal.hpp"
#include "test/predicate/istream_and_file_equal.hpp"

#include <fstream>
#include <sstream>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::boost::filesystem::current_path;
using ::boost::filesystem::path;
using ::std::ifstream;
using ::std::istream;
using ::std::istringstream;
using ::std::ostringstream;
using ::std::regex;
using ::std::string;
using ::std::stringstream;

namespace cath {
	namespace test {

		/// \brief The cath_superposer_test_suite_fixture to assist in testing cath_superposer
		struct cath_superposer_test_suite_fixture : protected global_test_constants {
		protected:
			~cath_superposer_test_suite_fixture() noexcept = default;

		public:
			const path MULTIPLE_B4DXN4_MODELS_FILE          { TEST_RESIDUE_IDS_DATA_DIR()   / "B4DXN4.multiple_END_separated_models.pdb" };
			const path CORRECT_B4DXN4_STDOUT_SUP_FILE       { TEST_RESIDUE_IDS_DATA_DIR()   / "B4DXN4.multiple_superposed_models.pdb"    };
			const path CORRECT_B4DXN4_SINGLE_GRAD_PYMOL_FILE{ TEST_RESIDUE_IDS_DATA_DIR()   / "B4DXN4.single_gradient_colour.pml"        };
			const path MULTIPLE_E9PB15_MODELS_FILE          { TEST_RESIDUE_IDS_DATA_DIR()   / "E9PB15.multiple_END_separated_models.pdb" };
			const path CORRECT_E9PB15_PYMOL_SUP_FILE        { TEST_RESIDUE_IDS_DATA_DIR()   / "E9PB15.multiple_superposed_models.pml"    };
			const path CORRECT_Q9HAU8_PYMOL_SUP_FILE        { TEST_RESIDUE_IDS_DATA_DIR()   / "Q9HAU8.multiple_superposed_models.pml"    };
			const path CORRECT_3_90_400_10_PYMOL_SUP_FILE   { TEST_MULTI_SSAP_SUPERPOSE_DIR() / "3.90.400.10.pml"                          };

			/// \brief The temporary output file
			const temp_file temp_cath_superposer_output_file{ ".temp_cath_superposer_test_file.%%%%-%%%%-%%%%-%%%%" };

			/// \brief The path of temporary output file
			const path temp_output_filename = get_filename( temp_cath_superposer_output_file );

			/// \brief The name of the cath_superpose executable
			///
			/// This doesn't really matter much and is just a dummy value in the tests
			const string CATH_SUPERPOSE_EXE = "cath-superpose";

			/// \brief An empty stdin stream
			istringstream empty_stdin;

			/// \brief One of the test directories
			const path orient_backbone_test_dir = TEST_SOURCE_DATA_DIR() / "orient" / "backbone_complete";

			void check_cath_superposer_use_case(const str_vec &,
			                                    istream &,
			                                    const path &,
			                                    const bool &);

			void check_cath_superposer_std_in_use_case(const str_vec &,
			                                           const path &,
			                                           const path &,
			                                           const bool &);
		};

	}  // namespace test
}  // namespace cath

void cath::test::cath_superposer_test_suite_fixture::check_cath_superposer_use_case(const str_vec &arg_command_line_args,    ///< TODOCUMENT
                                                                                    istream       &arg_istream,              ///< TODOCUMENT
                                                                                    const path    &arg_expected_output_file, ///< TODOCUMENT
                                                                                    const bool    &arg_outputs_to_temp_file  ///< TODOCUMENT
                                                                                    ) {
	argc_argv_faker faked_argc_and_argv(arg_command_line_args);

	const auto my_cath_superpose_options = make_and_parse_options<cath_superpose_options>(
		faked_argc_and_argv.get_argc(),
		faked_argc_and_argv.get_argv(),
		parse_sources::CMND_LINE_ONLY
	);

	// Prepare an ostringstream to capture the output
	stringstream test_stdout;
	stringstream test_stderr;

	// Perform the superposition and then close the input file
	cath_superposer::superpose(my_cath_superpose_options, arg_istream, test_stdout, test_stderr);

	if ( arg_outputs_to_temp_file ) {
		const auto output_file = temp_output_filename;
		if ( ! exists( output_file ) )  {
			BOOST_LOG_TRIVIAL( error ) << "cath-superpose command did not produce output file. Got stdout is: \""
			                           << test_stdout.str()
			                           << "\". Got stderr is: \""
			                           << test_stderr.str()
			                           << "\"";
		}

		// Blank out the version number in the superposition
		regex_replace_file(
			output_file,
			regex{ R"(cath-superpose \((v\d+\.\d+\.\d+\-\d+\-g?\w{8})?\))" },
			"cath-superpose (vX.X.X-X-XXXXXXXX)"
		);

		BOOST_CHECK_FILES_EQUAL(              output_file, arg_expected_output_file );
	}
	else {
		BOOST_CHECK_ISTREAM_AND_FILE_EQUAL(              test_stdout, "stdout_from_cath_superpose", arg_expected_output_file );
	}
}

void cath::test::cath_superposer_test_suite_fixture::check_cath_superposer_std_in_use_case(const str_vec &arg_command_line_args,    ///< TODOCUMENT
                                                                                           const path    &arg_input_file,           ///< TODOCUMENT
                                                                                           const path    &arg_expected_output_file, ///< TODOCUMENT
                                                                                           const bool    &arg_outputs_to_temp_file  ///< TODOCUMENT
                                                                                           ) {
	// Prepare an ifstream to read input from the input file
	ifstream stdin_like_ifstream;
	open_ifstream(stdin_like_ifstream, arg_input_file);
	check_cath_superposer_use_case(
		arg_command_line_args,
		stdin_like_ifstream,
		arg_expected_output_file,
		arg_outputs_to_temp_file
	);
	stdin_like_ifstream.close();
}

/// \TODO Remove hard-coded option name strings

BOOST_FIXTURE_TEST_SUITE(cath_superposer_test_suite, cath::test::cath_superposer_test_suite_fixture)

/// \brief
BOOST_AUTO_TEST_CASE(basic_genome3d_use_case) {
	check_cath_superposer_std_in_use_case(
		{ CATH_SUPERPOSE_EXE,
		  "--" + pdb_input_options_block::PO_PDBS_FROM_STDIN,
		  "--" + superposition_output_options_block::PO_SUP_TO_STDOUT,
		  "--" + alignment_input_options_block::PO_RES_NAME_ALIGN },
		MULTIPLE_B4DXN4_MODELS_FILE,
		CORRECT_B4DXN4_STDOUT_SUP_FILE,
		false
	);
}

/// \todo Make program exits work by throwing a special type of program_exit_exception that can be caught by
///       the test framework and then get the below test to check for exit with "No valid PDBs were loaded"
///// \brief
//BOOST_AUTO_TEST_CASE(genome3d_empty_stdin) {
//	check_cath_superposer_use_case(
//		{ CATH_SUPERPOSE_EXE,
//		  "--" + pdb_input_options_block::PO_PDBS_FROM_STDIN,
//		  "--" + superposition_output_options_block::PO_SUP_TO_STDOUT,
//		  "--" + alignment_input_options_block::PO_RES_NAME_ALIGN,
//		empty_stdin,
//		CORRECT_B4DXN4_STDOUT_SUP_FILE,
//		false
//	);
//}

/// \brief
BOOST_AUTO_TEST_CASE(genome3d_gradient_colour_alignment_use_case) {
	check_cath_superposer_std_in_use_case(
		{ CATH_SUPERPOSE_EXE,
		  "--" + pdb_input_options_block::PO_PDBS_FROM_STDIN,
		  "--" + superposition_output_options_block::PO_SUP_TO_PYMOL_FILE,
		  temp_output_filename.string(),
		  "--" + alignment_input_options_block::PO_RES_NAME_ALIGN,
		  "--gradient-colour-alignment" },
		MULTIPLE_E9PB15_MODELS_FILE,
		CORRECT_E9PB15_PYMOL_SUP_FILE,
		true
	);
}

/// \brief
BOOST_AUTO_TEST_CASE(genome3d_single_gradient_colour_alignment_use_case) {
	check_cath_superposer_std_in_use_case(
		{ CATH_SUPERPOSE_EXE,
		  "--" + pdb_input_options_block::PO_PDBS_FROM_STDIN,
		  "--" + superposition_output_options_block::PO_SUP_TO_PYMOL_FILE,
		  temp_output_filename.string(),
		  "--" + alignment_input_options_block::PO_RES_NAME_ALIGN,
		  "--gradient-colour-alignment" },
		TEST_RESIDUE_IDS_DATA_DIR() / "B4DXN4.DomSerf.1.pdb",
		CORRECT_B4DXN4_SINGLE_GRAD_PYMOL_FILE,
		true
	);
}

/// \brief
BOOST_AUTO_TEST_CASE(genome3d_big_tails_use_case) {
	check_cath_superposer_use_case(
		{ CATH_SUPERPOSE_EXE,
		  "--" + alignment_input_options_block::PO_RES_NAME_ALIGN,
		  "--pdb-infile",
		  ( TEST_RESIDUE_IDS_DATA_DIR() / "Q9HAU8.DomSerf.1.pdb"     ).string(),
		  "--pdb-infile",
		  ( TEST_RESIDUE_IDS_DATA_DIR() / "Q9HAU8.PHYRE2.1.pdb"      ).string(),
		  "--pdb-infile",
		  ( TEST_RESIDUE_IDS_DATA_DIR() / "Q9HAU8.SUPERFAMILY.2.pdb" ).string(),
		  "--pdb-infile",
		  ( TEST_RESIDUE_IDS_DATA_DIR() / "Q9HAU8.VIVACE.1.pdb"      ).string(),
		  "--" + superposition_output_options_block::PO_SUP_TO_PYMOL_FILE,
		  temp_output_filename.string(),
		  "--gradient-colour-alignment" },
		empty_stdin,
		CORRECT_Q9HAU8_PYMOL_SUP_FILE,
		true
	);
}

/// \brief
BOOST_AUTO_TEST_CASE(multi_ssap_3_90_400_10_use_case) {
	check_cath_superposer_use_case(
		{ CATH_SUPERPOSE_EXE,
		  "--ssap-scores-infile",
		  ( TEST_MULTI_SSAP_SUPERPOSE_DIR() / "3.90.400.10.ssap_scores_file" ).string(),
		  "--pdb-infile",
		  ( TEST_MULTI_SSAP_SUPERPOSE_DIR() / "1g5aA03" ).string(),
		  "--pdb-infile",
		  ( TEST_MULTI_SSAP_SUPERPOSE_DIR() / "1r7aA02" ).string(),
		  "--pdb-infile",
		  ( TEST_MULTI_SSAP_SUPERPOSE_DIR() / "1wzaA02" ).string(),
		  "--pdb-infile",
		  ( TEST_MULTI_SSAP_SUPERPOSE_DIR() / "1zjaA02" ).string(),
		  "--" + superposition_output_options_block::PO_SUP_TO_PYMOL_FILE,
		  temp_output_filename.string() },
		empty_stdin,
		CORRECT_3_90_400_10_PYMOL_SUP_FILE,
		true
	);
}

/// \brief Perform on testcase generated from subsets of PDBs 2a50, 3ako and 3cgl
///        that previously caused orient code to fail due to an attempt to access
///        the CA coord of a residue without a CA coord
///
/// \TODO Change cath-ssap's options to also use --pdb-infile (pdb_input_options_block::PO_PDB_INFILE)
///       so that explicit PDB files can be specified...
///       ...and then change the --do-the-ssaps code to use that to pass the explicit location of the file...
///       ...and then remove the current_path() stuff from this code...
///       ...and then remove the absolute() call in this test
BOOST_AUTO_TEST_CASE(can_orient_superposition_with_residues_that_are_not_backbone_complete) {
	// The expected data
	const path expected = absolute( orient_backbone_test_dir / "expected_supn.sup_pdb" );

	// Move to the directory containing the data (to make the SSAPs work correctly)
	const path orig_path = current_path();
	current_path( orient_backbone_test_dir );

	// Catch the logging output
	const stringstream_log_sink log_sink;
	try {
		check_cath_superposer_use_case(
			{
				CATH_SUPERPOSE_EXE,
				"--" + alignment_input_options_block::PO_DO_THE_SSAPS,
				"--" + pdb_input_options_block::PO_PDB_INFILE, "2a50",
				"--" + pdb_input_options_block::PO_PDB_INFILE, "3ako",
				"--" + pdb_input_options_block::PO_PDB_INFILE, "3cgl",
				"--" + superposition_output_options_block::PO_SUP_FILE, temp_output_filename.string()
			},
			empty_stdin,
			expected,
			true
		);
	}
	catch (...) {
		current_path( orig_path );
		throw;
	}
	current_path( orig_path );
}

BOOST_AUTO_TEST_SUITE_END()
