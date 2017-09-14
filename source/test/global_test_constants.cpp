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

#include "global_test_constants.hpp"

#include "exception/runtime_error_exception.hpp"
#include "ssap/ssap.hpp" // **** For reset_ssap_global_variables() ****

#include <limits>

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::filesystem::current_path;
using boost::filesystem::path;
using boost::filesystem::temp_directory_path;

constexpr double global_test_constants::DOUBLE_INFINITY;
constexpr double global_test_constants::DOUBLE_QUIET_NAN;
constexpr double global_test_constants::DOUBLE_SIGNALING_NAN;
constexpr double global_test_constants::DOUBLE_DENORM_MIN;
constexpr double global_test_constants::DOUBLE_MIN;
constexpr double global_test_constants::DOUBLE_MAX;

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
}

/// \brief Virtual empty destructor for global_test_constants.
global_test_constants::~global_test_constants() noexcept {
	try {
		cleanup_temporary_files();
	}
	catch (...) {
	}
}

// /// \brief Provide access to a static example residue_querier
// const residue_querier & global_test_constants::EXAMPLE_RESIDUE_QUERIER() {
// 	static const residue_querier example_residue_querier{};
// 	return example_residue_querier;
// }

// /// \brief Provide access to a static example sec_struc_querier
// const sec_struc_querier & global_test_constants::EXAMPLE_SEC_STRUC_QUERIER() {
// 	static const sec_struc_querier example_sec_struc_querier{};
// 	return example_sec_struc_querier;
// }

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
const path & global_test_constants::TEST_RESIDUE_IDS_DATA_DIR() {
	static const path test_residue_names_data_dir   ( TEST_SOURCE_DATA_DIR() / "residue_ids"   );
	return test_residue_names_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_SOURCE_DATA_DIR() {
	static const path test_source_data_dir          ( "build-test-data"             );
	return test_source_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_EXAMPLE_PDBS_DATA_DIR() {
	static const path test_example_pdbs_data_dir   ( TEST_SOURCE_DATA_DIR() / "example_pdbs" );
	return test_example_pdbs_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_SSAP_REGRESSION_DATA_DIR() {
	static const path test_ssap_regression_data_dir   ( TEST_SOURCE_DATA_DIR() / "ssap_regression" );
	return test_ssap_regression_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_SSAP_ALIGNMENT_GLUING_DATA_DIR() {
	static const path test_ssap_alignment_gluing_data_dir   ( TEST_SOURCE_DATA_DIR() / "ssap_scores_alignment_gluing" );
	return test_ssap_alignment_gluing_data_dir;
}

/// \brief Getter for directory of materials for testing superposition JSON
const path & global_test_constants::TEST_SUP_JSON_DIR() {
	static const path test_sup_json_data_dir{ TEST_SOURCE_DATA_DIR() / "superposition_json" };
	return test_sup_json_data_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_SVM_DIR() {
	static const path test_svm_dir{ TEST_SOURCE_DATA_DIR() / "svm" };
	return test_svm_dir;
}

/// \brief TODOCUMENT
const path & global_test_constants::TEST_OUTPUT_DIRECTORY() {
	static const path test_output_directory( temp_directory_path() );
	return test_output_directory;
}

/// \brief Test constant for the cath-resolve-hits test data directory
const path & global_test_constants::CRH_TEST_DATA_DIR() {
	static const path crh_test_data_dir( TEST_SOURCE_DATA_DIR() / "resolve_hits" );
	return crh_test_data_dir;
}

/// \brief Test constant for the cath-resolve-hits cath_gene3d_dc_handling test data subdirectory
const path & global_test_constants::CRH_CATH_DC_HANDLING_DATA_DIR() {
	static const path crh_cath_dc_handling_data_dir( CRH_TEST_DATA_DIR() / "cath_gene3d_dc_handling" );
	return crh_cath_dc_handling_data_dir;
}

/// \brief Test constant for the cath-resolve-hits hmm coverage test data subdirectory
const path & global_test_constants::CRH_HMM_COVERAGE_DATA_DIR() {
	static const path crh_hmm_coverage_data_dir( CRH_TEST_DATA_DIR() / "hmm_coverage" );
	return crh_hmm_coverage_data_dir;
}

/// \brief Test constant for the cath-resolve-hits hmmscan test data subdirectory
const path & global_test_constants::CRH_HMMSCAN_DATA_DIR() {
	static const path crh_hmmscan_data_dir( CRH_TEST_DATA_DIR() / "hmmscan" );
	return crh_hmmscan_data_dir;
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

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_DOMTBL_IN_FILENAME() {
	static const path crh_eg_domtbl_in_filename      ( CRH_TEST_DATA_DIR() / "eg_domtblout.in"                     );
	return crh_eg_domtbl_in_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_DOMTBL_OUT_FILENAME() {
	static const path crh_eg_domtbl_out_filename     ( CRH_TEST_DATA_DIR() / "eg_domtblout.out"                    );
	return crh_eg_domtbl_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_DOMTBL_LIMIT_2_OUT_FILENAME() {
	static const path crh_eg_domtbl_lim2_out_filename( CRH_TEST_DATA_DIR() / "eg_domtblout.limit_2.out"            );
	return crh_eg_domtbl_lim2_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_DOMTBL_JSON_OUT_FILENAME() {
	static const path crh_eg_domtbl_json_out_filename( CRH_TEST_DATA_DIR() / "eg_domtblout.json"                   );
	return crh_eg_domtbl_json_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_IN_FILENAME() {
	static const path crh_eg_hmmsearch_in_filename   ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.in"                 );
	return crh_eg_hmmsearch_in_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_out_filename              ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.out"                    );
	return crh_eg_hmmsearch_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_LIMIT_2_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_lim2_out_filename         ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.limit_2.out"            );
	return crh_eg_hmmsearch_lim2_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_BIG_GAP_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_big_gap_out_filename      ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.out_big_gap"            );
	return crh_eg_hmmsearch_big_gap_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_SMALL_GAP_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_small_gap_out_filename    ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.out_small_gap"          );
	return crh_eg_hmmsearch_small_gap_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_TRIMMED_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_trimmed_out_filename      ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.out_output_trimmed"     );
	return crh_eg_hmmsearch_trimmed_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_BIG_TRIM_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_big_trim_out_filename     ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.out_big_trim"           );
	return crh_eg_hmmsearch_big_trim_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_HMMSEARCH_ALN_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_hmmsearch_aln_out_filename( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.out_hmmsearch_aln_summ" );
	return crh_eg_hmmsearch_hmmsearch_aln_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_SUMMARISE_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_summarise_filename        ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.summarise_out"          );
	return crh_eg_hmmsearch_summarise_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_HTML_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_html_filename             ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.html"                   );
	return crh_eg_hmmsearch_html_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_HMMSEARCH_JSON_OUT_FILENAME() {
	static const path crh_eg_hmmsearch_json_filename             ( CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.json"                   );
	return crh_eg_hmmsearch_json_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_RAW_EVALUE_IN_FILENAME() {
	static const path crh_eg_raw_eva_in_filename     ( CRH_TEST_DATA_DIR() / "eg_raw_evalue.in"                    );
	return crh_eg_raw_eva_in_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_RAW_EVALUE_OUT_FILENAME() {
	static const path crh_eg_raw_eva_out_filename    ( CRH_TEST_DATA_DIR() / "eg_raw_evalue.out"                   );
	return crh_eg_raw_eva_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_RAW_EVALUE_LIMIT_2_OUT_FILENAME() {
	static const path crh_eg_raw_eva_lim2_out_filename( CRH_TEST_DATA_DIR() / "eg_raw_evalue.limit_2.out"          );
	return crh_eg_raw_eva_lim2_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_RAW_SCORE_IN_FILENAME() {
	static const path crh_eg_raw_sc_in_filename      ( CRH_TEST_DATA_DIR() / "eg_raw_score.in"                     );
	return crh_eg_raw_sc_in_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_RAW_SCORE_OUT_FILENAME() {
	static const path crh_eg_raw_sc_out_filename     ( CRH_TEST_DATA_DIR() / "eg_raw_score.out"                    );
	return crh_eg_raw_sc_out_filename;
}

/// \brief Test constant for a cath-resolve-hits test file
const path & global_test_constants::CRH_EG_RAW_SCORE_LIMIT_2_OUT_FILENAME() {
	static const path crh_eg_raw_sc_lim2_out_filename( CRH_TEST_DATA_DIR() / "eg_raw_score.limit_2.out"            );
	return crh_eg_raw_sc_lim2_out_filename;
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
	static constexpr float looser_accuracy_percentage_float = static_cast<float>( 1.0 );
	return looser_accuracy_percentage_float;
}

/// \brief TODOCUMENT
template <>
const double & global_test_constants::LOOSER_ACCURACY_PERCENTAGE_TMPL<double>() {
	static constexpr double looser_accuracy_percentage_double ( 0.0001 );
	return looser_accuracy_percentage_double;
}

/// \brief TODOCUMENT
template <>
const long double & global_test_constants::LOOSER_ACCURACY_PERCENTAGE_TMPL<long double>() {
	static constexpr long double looser_accuracy_percentage_long_double ( 0.00001 );
	return looser_accuracy_percentage_long_double;
}

/// \brief TODOCUMENT
template <>
const float & global_test_constants::ACCURACY_PERCENTAGE_TMPL<float>() {
	static constexpr float accuracy_percentage = static_cast<float>( 0.0000001 );
	return accuracy_percentage;
}

/// \brief TODOCUMENT
template <>
const double & global_test_constants::ACCURACY_PERCENTAGE_TMPL<double>() {
	static constexpr double accuracy_percentage  ( 0.000000001                             );
	return accuracy_percentage;
}

/// \brief TODOCUMENT
template <>
const long double & global_test_constants::ACCURACY_PERCENTAGE_TMPL<long double>() {
	static constexpr long double accuracy_percentage  ( 0.00000000001                             );
	return accuracy_percentage;
}

/// \brief TODOCUMENT
const doub_vec & global_test_constants::INVALID_DOUBLES() {
	static const doub_vec invalid_doubles = {
		DOUBLE_INFINITY,
		DOUBLE_QUIET_NAN,
		DOUBLE_SIGNALING_NAN
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

