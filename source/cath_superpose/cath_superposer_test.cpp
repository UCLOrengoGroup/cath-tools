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

#include <boost/log/trivial.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "cath_superposer.h"
#include "common/argc_argv_faker.h"
#include "common/file/open_fstream.h"
#include "common/test_predicate/files_equal.h"
#include "common/test_predicate/istream_and_file_equal.h"
#include "common/file/temp_file.h"
#include "options/executable/cath_superpose_options/cath_superpose_options.h"
#include "test/global_test_constants.h"

#include <fstream>
#include <sstream>

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The cath_superposer_test_suite_fixture to assist in testing cath_superposer
		struct cath_superposer_test_suite_fixture : protected global_test_constants {
		protected:
			~cath_superposer_test_suite_fixture() noexcept = default;

		public:
			const path MULTIPLE_B4DXN4_MODELS_FILE          { TEST_RESIDUE_NAMES_DATA_DIR()   / "B4DXN4.multiple_END_separated_models.pdb" };
			const path CORRECT_B4DXN4_STDOUT_SUP_FILE       { TEST_RESIDUE_NAMES_DATA_DIR()   / "B4DXN4.multiple_superposed_models.pdb"    };
			const path CORRECT_B4DXN4_SINGLE_GRAD_PYMOL_FILE{ TEST_RESIDUE_NAMES_DATA_DIR()   / "B4DXN4.single_gradient_colour.pml"        };
			const path MULTIPLE_E9PB15_MODELS_FILE          { TEST_RESIDUE_NAMES_DATA_DIR()   / "E9PB15.multiple_END_separated_models.pdb" };
			const path CORRECT_E9PB15_PYMOL_SUP_FILE        { TEST_RESIDUE_NAMES_DATA_DIR()   / "E9PB15.multiple_superposed_models.pml"    };
			const path CORRECT_Q9HAU8_PYMOL_SUP_FILE        { TEST_RESIDUE_NAMES_DATA_DIR()   / "Q9HAU8.multiple_superposed_models.pml"    };
			const path CORRECT_3_90_400_10_PYMOL_SUP_FILE   { TEST_MULTI_SSAP_SUPERPOSE_DIR() / "3.90.400.10.pml"                          };

			const temp_file temp_cath_superposer_output_file{ ".temp_cath_superposer_test_file.%%%%-%%%%-%%%%-%%%%" };

			void check_cath_superposer_use_case(const str_vec &,
			                                    istream &,
			                                    const path &,
			                                    const bool &);

			void check_cath_superposer_std_in_use_case(const str_vec &,
			                                           const path &,
			                                           const path &,
			                                           const bool &);
		};

	}
}

void cath::test::cath_superposer_test_suite_fixture::check_cath_superposer_use_case(const str_vec &arg_command_line_args,    ///< TODOCUMENT
                                                                                    istream       &arg_istream,              ///< TODOCUMENT
                                                                                    const path    &arg_expected_output_file, ///< TODOCUMENT
                                                                                    const bool    &arg_outputs_to_temp_file  ///< TODOCUMENT
                                                                                    ) {
	argc_argv_faker faked_argc_and_argv(arg_command_line_args);
	cath_superpose_options my_cath_superpose_options;
	my_cath_superpose_options.parse_options(
		faked_argc_and_argv.get_argc(),
		faked_argc_and_argv.get_argv()
	);

	// Prepare an ostringstream to capture the output
	stringstream test_stdout;
	stringstream test_stderr;

	// Perform the superposition and then close the input file
	cath_superposer::superpose(my_cath_superpose_options, arg_istream, test_stdout, test_stderr);

	if ( arg_outputs_to_temp_file ) {
		const auto output_file = get_filename( temp_cath_superposer_output_file );
		if ( ! exists( output_file ) )  {
			BOOST_LOG_TRIVIAL( error ) << "cath-superpose command did not produce output file. Got stdout is: \""
			                           << test_stdout.str()
									   << "\". Got stderr is: \""
									   << test_stderr.str()
									   << "\"";
		}
		BOOST_CHECK_FILES_EQUAL(              output_file, arg_expected_output_file );
//		BOOST_CHECK_FILES_EQUAL_OR_OVERWRITE( output_file, arg_expected_output_file );
	}
	else {
		BOOST_CHECK_ISTREAM_AND_FILE_EQUAL(              test_stdout, "stdout_from_cath_superpose", arg_expected_output_file );
//		BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( test_stdout, "stdout_from_cath_superpose", arg_expected_output_file );
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


BOOST_FIXTURE_TEST_SUITE(cath_superposer_test_suite, cath::test::cath_superposer_test_suite_fixture)

/// \brief
BOOST_AUTO_TEST_CASE(basic_genome3d_use_case) {
	check_cath_superposer_std_in_use_case(
		{ "cath-superpose",
		  "--pdbs-from-stdin",
		  "--sup-to-stdout",
		  "--res-name-align" },
		MULTIPLE_B4DXN4_MODELS_FILE,
		CORRECT_B4DXN4_STDOUT_SUP_FILE,
		false
	);
}

/// \todo Make program exits work by throwing a special type of program_exit_exception that can be caught by
///       the test framework and then get the below test to check for exit with "No valid PDBs were loaded"
///// \brief
//BOOST_AUTO_TEST_CASE(genome3d_empty_stdin) {
//	istringstream empty_stdin;
//	check_cath_superposer_use_case(
//		{ "cath-superpose",
//		  "--pdbs-from-stdin",
//		  "--sup-to-stdout",
//		  "--res-name-align" },
//		empty_stdin,
//		CORRECT_B4DXN4_STDOUT_SUP_FILE,
//		false
//	);
//}

/// \brief
BOOST_AUTO_TEST_CASE(genome3d_gradient_colour_alignment_use_case) {
	check_cath_superposer_std_in_use_case(
		{ "cath-superpose",
		  "--pdbs-from-stdin",
		  "--sup-to-pymol-file",
		  get_filename( temp_cath_superposer_output_file ).string(),
		  "--res-name-align",
		  "--gradient-colour-alignment" },
		MULTIPLE_E9PB15_MODELS_FILE,
		CORRECT_E9PB15_PYMOL_SUP_FILE,
		true
	);
}

/// \brief
BOOST_AUTO_TEST_CASE(genome3d_single_gradient_colour_alignment_use_case) {
	check_cath_superposer_std_in_use_case(
		{ "cath-superpose",
		  "--pdbs-from-stdin",
		  "--sup-to-pymol-file",
		  get_filename( temp_cath_superposer_output_file ).string(),
		  "--res-name-align",
		  "--gradient-colour-alignment" },
		TEST_RESIDUE_NAMES_DATA_DIR() / "B4DXN4.DomSerf.1.pdb",
		CORRECT_B4DXN4_SINGLE_GRAD_PYMOL_FILE,
		true
	);
}

/// \brief
BOOST_AUTO_TEST_CASE(genome3d_big_tails_use_case) {
	istringstream empty_stdin;
	check_cath_superposer_use_case(
		{ "cath-superpose",
		  "--res-name-align",
		  "--pdb-infile",
		  ( TEST_RESIDUE_NAMES_DATA_DIR() / "Q9HAU8.DomSerf.1.pdb"     ).string(),
		  "--pdb-infile",
		  ( TEST_RESIDUE_NAMES_DATA_DIR() / "Q9HAU8.PHYRE2.1.pdb"      ).string(),
		  "--pdb-infile",
		  ( TEST_RESIDUE_NAMES_DATA_DIR() / "Q9HAU8.SUPERFAMILY.2.pdb" ).string(),
		  "--pdb-infile",
		  ( TEST_RESIDUE_NAMES_DATA_DIR() / "Q9HAU8.VIVACE.1.pdb"      ).string(),
		  "--sup-to-pymol-file",
		  get_filename( temp_cath_superposer_output_file ).string(),
		  "--gradient-colour-alignment" },
		empty_stdin,
		CORRECT_Q9HAU8_PYMOL_SUP_FILE,
		true
	);
}

/// \brief
BOOST_AUTO_TEST_CASE(multi_ssap_3_90_400_10_use_case) {
	istringstream empty_stdin;
	check_cath_superposer_use_case(
		{ "cath-superpose",
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
		  "--sup-to-pymol-file",
		  get_filename( temp_cath_superposer_output_file ).string() },
		empty_stdin,
		CORRECT_3_90_400_10_PYMOL_SUP_FILE,
		true
	);
}

BOOST_AUTO_TEST_SUITE_END()
