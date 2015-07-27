/// \file
/// \brief The global_test_constants class definitions

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

#include "global_test_constants.h"

#include "exception/runtime_error_exception.h"
#include "ssap/ssap.h" // **** For reset_ssap_global_variables() ****

//#include <iostream> // *** TEMPORARY ***
#include <limits>

using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace std;

/// \brief TODOCUMENT
void global_test_constants::check_required_files_exist() {
	for (const path &required_file : REQUIRED_PREEXISTING_FILES() ) {
		if ( ! exists( required_file ) ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				"Unable to find required test file \""
				+ required_file.string()
				+ "\" (current directory is \""
				+ current_path().string()
				+ "\")"
			));
		}
	}
}

/// \brief TODOCUMENT
void global_test_constants::cleanup_temporary_files() {
	for (const path &temporary_file : TEMPORARY_FILES()) {
		if ( exists( temporary_file ) ) {
			remove( temporary_file );
		}
		if ( exists( NONEXISTENT_FILE() ) ) {
			remove( NONEXISTENT_FILE() );
		}
	}
}

/// \brief Default ctor for global_test_constants.
global_test_constants::global_test_constants() {
	check_required_files_exist();
	cleanup_temporary_files();
	reset_ssap_global_variables();
}

/// \brief Virtual empty destructor for global_test_constants.
global_test_constants::~global_test_constants() noexcept {
	try {
		cleanup_temporary_files();
	}
	catch (...) {
	}
}

/// \brief Provide access to a static example residue_querier
const residue_querier & global_test_constants::EXAMPLE_RESIDUE_QUERIER() {
	static const residue_querier example_residue_querier{};
	return example_residue_querier;
}

/// \brief Provide access to a static example sec_struc_querier
const sec_struc_querier & global_test_constants::EXAMPLE_SEC_STRUC_QUERIER() {
	static const sec_struc_querier example_sec_struc_querier{};
	return example_sec_struc_querier;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_BASIC_FILE_TEST_DATA_DIR() {
	static const path test_basic_file_test_data_dir ( TEST_SOURCE_DATA_DIR() / "basic_file_test" );
	return test_basic_file_test_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_MULTI_SSAP_SUPERPOSE_DIR() {
	static const path test_basic_file_test_data_dir ( TEST_SOURCE_DATA_DIR() / "multi_ssap_superpose_3.90.400.10" );
	return test_basic_file_test_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_RESIDUE_NAMES_DATA_DIR() {
	static const path test_residue_names_data_dir   ( TEST_SOURCE_DATA_DIR() / "residue_names"   );
	return test_residue_names_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_SOURCE_DATA_DIR() {
	static const path test_source_data_dir          ( "build-test-data"             );
	return test_source_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_SSAP_REGRESSION_DATA_DIR() {
	static const path test_ssap_regression_data_dir   ( TEST_SOURCE_DATA_DIR() / "ssap_regression" );
	return test_ssap_regression_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_OUTPUT_DIRECTORY() {
	static const path test_output_directory( temp_directory_path() );
	return test_output_directory;
}

/// \brief TODOCUMENT
const path & global_test_constants::NONEXISTENT_FILE() {
	static const path nonexistent_file( "filename_of_nonexistent_file" );
	return nonexistent_file;
}

/// \brief TODOCUMENT
const string & global_test_constants::EXAMPLE_A_PDB_STEMNAME() {
	static const string example_a_pdb_stemname( "1c0pA01" );
	return example_a_pdb_stemname;
}

/// \brief TODOCUMENT
const string & global_test_constants::EXAMPLE_B_PDB_STEMNAME() {
	static const string example_b_pdb_stemname( "1hdoA00" );
	return example_b_pdb_stemname;
}

/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_A_PDB_FILENAME() {
	static const path example_a_pdb_filename  ( TEST_SOURCE_DATA_DIR()  / EXAMPLE_A_PDB_STEMNAME()               );
	return example_a_pdb_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_B_PDB_FILENAME() {
	static const path example_b_pdb_filename  ( TEST_SOURCE_DATA_DIR()  / EXAMPLE_B_PDB_STEMNAME()               );
	return example_b_pdb_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_A_DSSP_FILENAME() {
	static const path example_a_dssp_filename  ( TEST_SOURCE_DATA_DIR()  / ( EXAMPLE_A_PDB_STEMNAME() + ".dssp" ) );
	return example_a_dssp_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_B_DSSP_FILENAME() {
	static const path example_b_dssp_filename  ( TEST_SOURCE_DATA_DIR()  / ( EXAMPLE_B_PDB_STEMNAME() + ".dssp" ) );
	return example_b_dssp_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_A_WOLF_FILENAME() {
	static const path example_a_wolf_filename ( TEST_SOURCE_DATA_DIR()  / ( EXAMPLE_A_PDB_STEMNAME() + ".wolf" ) );
	return example_a_wolf_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_B_WOLF_FILENAME() {
	static const path example_b_wolf_filename ( TEST_SOURCE_DATA_DIR()  / ( EXAMPLE_B_PDB_STEMNAME() + ".wolf" ) );
	return example_b_wolf_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_A_SEC_FILENAME() {
	static const path example_a_sec_filename  ( TEST_SOURCE_DATA_DIR()  / ( EXAMPLE_A_PDB_STEMNAME() + ".sec"  ) );
	return example_a_sec_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_B_SEC_FILENAME() {
	static const path example_b_sec_filename  ( TEST_SOURCE_DATA_DIR()  / ( EXAMPLE_B_PDB_STEMNAME() + ".sec"  ) );
	return example_b_sec_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::ALIGNMENT_FILE() {
	static const path alignment_file        ( TEST_SOURCE_DATA_DIR()  / "1c0pA011hdoA00.list"               );
	return alignment_file;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_PAIR_SUPPOSN_XML() {
	static const path test_pair_supposn_xml ( TEST_SOURCE_DATA_DIR()  / "simple_pairwise_superposition.xml" );
	return test_pair_supposn_xml;
}

/// \brief TODOCUMENT
const path & global_test_constants::MODIFIED_PDB_FILENAME() {
	static const path modified_pdb_filename ( TEST_OUTPUT_DIRECTORY() / "1c0pA01.modified"                  );
	return modified_pdb_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_OUTPUT_FILENAME() {
	static const path test_output_filename  ( TEST_OUTPUT_DIRECTORY() / "temp_test_1c0pA01.modified"        );
	return test_output_filename;
}

/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_DOUBLES_FILENAME() {
	static const path test_output_filename  ( TEST_SOURCE_DATA_DIR() / "simple_file_read_write" / "example_numbers" );
	return test_output_filename;
}
/// \brief TODOCUMENT
const path & global_test_constants::EXAMPLE_TUPLES_FILENAME() {
	static const path test_output_filename  ( TEST_SOURCE_DATA_DIR() / "simple_file_read_write" / "example_tuples" );
	return test_output_filename;
}

/// \brief TODOCUMENT
const double & global_test_constants::LOOSER_ACCURACY_PERCENTAGE() {
	return LOOSER_ACCURACY_PERCENTAGE_TMPL<double>();
}

/// \brief TODOCUMENT
const double & global_test_constants::ACCURACY_PERCENTAGE() {
	return ACCURACY_PERCENTAGE_TMPL<double>();
}

/// \brief TODOCUMENT
template <>
const float & global_test_constants::LOOSER_ACCURACY_PERCENTAGE_TMPL<float>() {
	static const float looser_accuracy_percentage_float = static_cast<float>( 1.0 );
	return looser_accuracy_percentage_float;
}

/// \brief TODOCUMENT
template <>
const double & global_test_constants::LOOSER_ACCURACY_PERCENTAGE_TMPL<double>() {
	static const double looser_accuracy_percentage_double ( 0.0001 );
	return looser_accuracy_percentage_double;
}

/// \brief TODOCUMENT
template <>
const long double & global_test_constants::LOOSER_ACCURACY_PERCENTAGE_TMPL<long double>() {
	static const long double looser_accuracy_percentage_long_double ( 0.00001 );
	return looser_accuracy_percentage_long_double;
}

/// \brief TODOCUMENT
template <>
const float & global_test_constants::ACCURACY_PERCENTAGE_TMPL<float>() {
	static const float accuracy_percentage = static_cast<float>( 0.0000001 );
	return accuracy_percentage;
}

/// \brief TODOCUMENT
template <>
const double & global_test_constants::ACCURACY_PERCENTAGE_TMPL<double>() {
	static const double accuracy_percentage  ( 0.000000001                             );
	return accuracy_percentage;
}

/// \brief TODOCUMENT
template <>
const long double & global_test_constants::ACCURACY_PERCENTAGE_TMPL<long double>() {
	static const long double accuracy_percentage  ( 0.00000000001                             );
	return accuracy_percentage;
}

/// \brief TODOCUMENT
const double & global_test_constants::DOUBLE_INFINITY() {
	static const double double_infinity      ( numeric_limits<double>::infinity()      );
	return double_infinity;
}

/// \brief TODOCUMENT
const double & global_test_constants::DOUBLE_QUIET_NAN() {
	static const double double_quiet_nan     ( numeric_limits<double>::quiet_NaN()     );
	return double_quiet_nan;
}

/// \brief TODOCUMENT
const double & global_test_constants::DOUBLE_SIGNALING_NAN() {
	static const double double_signaling_nan ( numeric_limits<double>::signaling_NaN() );
	return double_signaling_nan;
}

/// \brief TODOCUMENT
const double & global_test_constants::DOUBLE_DENORM_MIN() {
	static const double double_denorm_min    ( numeric_limits<double>::denorm_min()    );
	return double_denorm_min;
}

/// \brief TODOCUMENT
const double & global_test_constants::DOUBLE_MIN() {
	static const double double_min           ( numeric_limits<double>::min()           );
	return double_min;
}

/// \brief TODOCUMENT
const double & global_test_constants::DOUBLE_MAX() {
	static const double double_max           ( numeric_limits<double>::max()           );
	return double_max;
}
/// \brief TODOCUMENT
const doub_vec & global_test_constants::INVALID_DOUBLES() {
	static const doub_vec invalid_doubles = {
		DOUBLE_INFINITY(),
		DOUBLE_QUIET_NAN(),
		DOUBLE_SIGNALING_NAN()
	};
	return invalid_doubles;
}

/// \brief TODOCUMENT
const path_vec & global_test_constants::REQUIRED_PREEXISTING_FILES() {
	static const path_vec required_preexisting_files = {
		EXAMPLE_A_PDB_FILENAME(),
		EXAMPLE_B_PDB_FILENAME(),
		ALIGNMENT_FILE(),
		EXAMPLE_DOUBLES_FILENAME(),
		EXAMPLE_TUPLES_FILENAME()
	};
	return required_preexisting_files;
}

/// \brief TODOCUMENT
const path_vec & global_test_constants::TEMPORARY_FILES() {
	static const path_vec temporary_files = {
		MODIFIED_PDB_FILENAME(),
		TEST_OUTPUT_FILENAME()
	};
	return temporary_files;
}

