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

#include <boost/test/auto_unit_test.hpp>

#include "cath_refine_align/options/cath_refine_align_options.h"
#include "cath_superpose/cath_superposer.h"
#include "common/argc_argv_faker.h"
#include "common/file/open_fstream.h"
#include "common/file/temp_file.h"
#include "common/test_predicate/files_equal.h"
#include "common/test_predicate/istream_and_file_equal.h"
#include "test/global_test_constants.h"

#include <fstream>
#include <sstream>

using namespace boost::filesystem;
using namespace cath::common;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The cath_align_refiner_test_suite_fixture to assist in testing cath_align_refiner
		struct cath_align_refiner_test_suite_fixture : protected global_test_constants {
		protected:
			~cath_align_refiner_test_suite_fixture() noexcept = default;

			const path B4DXN4_MULTIPLE_MODELS_FILE    { TEST_RESIDUE_NAMES_DATA_DIR() / "B4DXN4.multiple_END_separated_models.pdb" };
			const path B4DXN4_CORRECT_STDOUT_SUP_FILE { TEST_RESIDUE_NAMES_DATA_DIR() / "B4DXN4.multiple_superposed_models.pdb"    };
			const path E9PB15_MULTIPLE_MODELS_FILE    { TEST_RESIDUE_NAMES_DATA_DIR() / "E9PB15.multiple_END_separated_models.pdb" };
			const path E9PB15_CORRECT_PYMOL_SUP_FILE  { TEST_RESIDUE_NAMES_DATA_DIR() / "E9PB15.multiple_superposed_models.pml"    };

			const temp_file temp_cath_align_refiner_output_file{ ".temp_cath_align_refiner_test_file.%%%%-%%%%-%%%%-%%%%" };

			void check_cath_align_refiner_use_case(const str_vec &,
			                                       istream &,
			                                       const path &,
			                                       const bool &);

			void check_cath_align_refiner_std_in_use_case(const str_vec &,
			                                              const path &,
			                                              const path &,
			                                              const bool &);
		};

	}
}

void cath::test::cath_align_refiner_test_suite_fixture::check_cath_align_refiner_use_case(const str_vec &arg_command_line_args,    ///< TODOCUMENT
                                                                                          istream       &arg_istream,              ///< TODOCUMENT
                                                                                          const path    &arg_expected_output_file, ///< TODOCUMENT
                                                                                          const bool    &arg_outputs_to_temp_file  ///< TODOCUMENT
                                                                                          ) {
		argc_argv_faker faked_argc_and_argv( arg_command_line_args );
		const auto my_cath_refine_align_options = make_and_parse_options<cath_refine_align_options>(
			faked_argc_and_argv.get_argc(),
			faked_argc_and_argv.get_argv()
		);

		// Prepare an ostringstream to capture the output
		stringstream test_ostream;

		// Perform the superposition and then close the input file
		cath_superposer::superpose( my_cath_refine_align_options, arg_istream, test_ostream );

		if (arg_outputs_to_temp_file) {
			BOOST_CHECK_FILES_EQUAL( get_filename( temp_cath_align_refiner_output_file ), arg_expected_output_file);
		}
		else {
			BOOST_CHECK_ISTREAM_AND_FILE_EQUAL(test_ostream, "stdout_from_cath_align_refine", arg_expected_output_file);
		}
	}

void cath::test::cath_align_refiner_test_suite_fixture::check_cath_align_refiner_std_in_use_case(const str_vec &arg_command_line_args,    ///< TODOCUMENT
                                                                                                 const path    &arg_input_file,           ///< TODOCUMENT
                                                                                                 const path    &arg_expected_output_file, ///< TODOCUMENT
                                                                                                 const bool    &arg_outputs_to_temp_file  ///< TODOCUMENT
                                                                                                 ) {
	// Prepare an ifstream to read input from the input file
	ifstream stdin_like_ifstream;
	open_ifstream(stdin_like_ifstream, arg_input_file);
	check_cath_align_refiner_use_case(
		arg_command_line_args,
		stdin_like_ifstream,
		arg_expected_output_file,
		arg_outputs_to_temp_file
	);
	stdin_like_ifstream.close();
};

BOOST_FIXTURE_TEST_SUITE(cath_align_refiner_test_suite, cath::test::cath_align_refiner_test_suite_fixture)

/// \brief
BOOST_AUTO_TEST_CASE(basic_genome3d_use_case) {
	check_cath_align_refiner_std_in_use_case(
		{
			"cath_align_refine",
			"--pdbs-from-stdin",
			"--sup-to-stdout",
			"--res-name-align"
		},
		B4DXN4_MULTIPLE_MODELS_FILE,
		B4DXN4_CORRECT_STDOUT_SUP_FILE,
		false
	);
}

/// \brief
BOOST_AUTO_TEST_CASE(genome3d_gradient_colour_alignment_use_case) {
	check_cath_align_refiner_std_in_use_case(
		{
			"cath_align_refine",
			"--pdbs-from-stdin",
			"--sup-to-pymol-file",
			get_filename( temp_cath_align_refiner_output_file ).string(),
			"--res-name-align",
			"--gradient-colour-alignment"
		},
		E9PB15_MULTIPLE_MODELS_FILE,
		E9PB15_CORRECT_PYMOL_SUP_FILE,
		true
	);
}

BOOST_AUTO_TEST_SUITE_END()
