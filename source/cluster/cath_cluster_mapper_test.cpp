/// \file
/// \brief The cath_cluster_mapper test suite

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

#include <boost/range/join.hpp>
#include <boost/test/auto_unit_test.hpp>
// #include <boost/algorithm/string/predicate.hpp>

#include "cluster/cath_cluster_mapper.hpp"
#include "cluster/options/options_block/clust_mapping_options_block.hpp"
#include "cluster/options/options_block/clustmap_input_options_block.hpp"
#include "cluster/options/options_block/clustmap_output_options_block.hpp"
#include "cluster/test/map_clusters_fixture.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/file/temp_file.hpp"
#include "common/test_predicate/files_equal.hpp"
#include "common/type_aliases.hpp"
// #include "cluster/options/options_block/clust_mapping_options_block.hpp"
// #include "cluster/options/options_block/clust_mapping_options_block.hpp"
// #include "cluster/options/options_block/clust_mapping_options_block.hpp"
// #include "cluster/options/options_block/clust_mapping_options_block.hpp"
// #include "cluster/options/options_block/clust_mapping_options_block.hpp"
// #include "cluster/options/options_block/clustmap_input_options_block.hpp"
// #include "cluster/options/options_block/clustmap_output_options_block.hpp"
// #include "cluster/options/options_block/crh_score_options_block.hpp"
// #include "cluster/test/map_clusters_fixture.hpp"
// #include "common/boost_addenda/log/log_to_ostream_guard.hpp"
// #include "common/boost_addenda/log/log_to_ostream_guard.hpp"
// #include "common/boost_addenda/test/boost_check_no_throw_diag.hpp"
// #include "common/file/read_string_from_file.hpp"
// #include "common/file/simple_file_read_write.hpp"
// #include "common/file/temp_file.hpp"
// #include "common/test_predicate/istream_and_file_equal.hpp"
// #include "test/global_test_constants.hpp"

#include <regex>
#include <sstream>

namespace cath { namespace test { } }

// using namespace cath::rslv;
// using namespace cath;
using namespace cath::clust;
using namespace cath::common;
using namespace cath::test;
using namespace std::literals::string_literals;

// using boost::algorithm::contains;
// using cath::common::copy_build;
// using cath::common::temp_file;
// using cath::common::write_file;
// using std::string;
using boost::filesystem::path;
using boost::range::join;
using std::istringstream;
using std::regex;
using std::ostringstream;

namespace cath {
	namespace test {

		/// \brief The cluster_mapper_test_suite_fixture to assist in testing calc_hit_list
		struct cluster_mapper_test_suite_fixture : protected map_clusters_fixture {
		protected:
			~cluster_mapper_test_suite_fixture() noexcept  = default;

			/// \brief Call perform_map_clusters() with the specified arguments preceded by a pseudo-program-name
			///        and the fixture's i/o streams.
			void execute_perform_map_clusters(const str_vec &arg_arguments ///< The arguments to pass to perform_map_clusters(), preceded by a pseudo-program-name
			                                  ) {
				const auto progname_rng = { "pseudo_program_name"s };
				perform_map_clusters(
					copy_build<str_vec>( join(
						progname_rng,
						arg_arguments
					) ),
					input_ss, output_ss
				);
			}

			/// \brief The input stream to use in the tests
			istringstream   input_ss;

			/// \brief The output stream to use in the tests
			ostringstream   output_ss;

			// /// \brief An output stream to which logging can be sent
			// ostringstream   log_ss;

			/// \brief A temporary temp_file
			const temp_file TEMP_TEST_FILE{ ".cath_cluster_mapper__test_file.%%%%-%%%%-%%%%-%%%%" };

			/// \brief The path of the TEMP_TEST_FILE
			const path      TEMP_TEST_FILE_FILENAME = get_filename( TEMP_TEST_FILE );
		};

	}  // namespace test
}  // namespace cath

/// \TODO: Tests to add:
///  * Gives helpful error on both --map-from-clustmemb-file and --read-batches-from-input specified
///  * Allows any string for input data cluster names
///  * Rejects --map-from-clustmemb-file cluster names that aren't positive-integers
///  * Overlapping domains are rejected
///  * Exact duplicates are warned *once*.
///  * Duplicates spotted if specifying in different ways (ie /1-2,3-4 vs 1-2_3-4)
///  * Accepts either overlap value being 50 or 100
///  * Rejects overlap value < 50 or > 100
///  * Enforce that both overlap values must be ≥ 50 and ≤ 100 but can be either 50 or 100.
///  * Errors if no input option specified
///  * Errors on attempt to specify --summary when not performing mappings

BOOST_FIXTURE_TEST_SUITE(cath_cluster_mapper_test_suite, cluster_mapper_test_suite_fixture)

// BOOST_AUTO_TEST_CASE(fails_on_attempt_to_mix_deprecated_options_with_new) {
// 	// When calling perform_map_clusters with options: - (a dash to read from the input stream)
// 	execute_perform_map_clusters( {
// 		"--" + clustmap_output_options_block::PO_JSON_OUTPUT_TO_FILE, "-"
// 		"--" + clust_mapping_options_block::PO_SUMMARISE,
// 		"--" + clust_mapping_options_block::PO_OUTPUT_FILE, TEMP_TEST_FILE_FILENAME.string(),
// 		"--" + clustmap_output_options_block::PO_QUIET,
// 	} );

// 	// Then expect the correct error message in the output stream
// 	BOOST_CHECK_EQUAL( output_ss.str(), "cath-resolve-hits: Cannot mix old, deprecated options (--output-file) with new, replacement options (--quiet, --json-output-to-file). Please use the new options only.\nSee 'cath-resolve-hits --help' for usage.\n" );
// }


BOOST_AUTO_TEST_CASE(fails_if_given_no_options) {
	// Given an input stream containing the input
	input_ss.str( eg_input_str() );

	// When calling perform_map_clusters with options: - (a dash to read from the input stream)
	execute_perform_map_clusters( { } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(Must specify an input file)" } ) );
}


BOOST_AUTO_TEST_CASE(processes_from_stdin_to_stdout) {
	// Given an input stream containing the input
	input_ss.str( eg_input_str() );

	// When calling perform_map_clusters with options: - (a dash to read from the input stream)
	execute_perform_map_clusters( { "-" } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_EQUAL( output_ss.str(), eg_renumber_only_result_str() );
}

BOOST_AUTO_TEST_CASE(processes_from_file_to_stdout) {
	// When calling perform_map_clusters with options: an input file
	execute_perform_map_clusters( { eg_input_file().string() } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_EQUAL( output_ss.str(), eg_renumber_only_result_str() );
}

// clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE
// clustmap_input_options_block::PO_READ_BATCHES_FROM_INPUT
//
// clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL
// clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL
//
// clustmap_input_options_block::PO_APPEND_BATCH_ID
// clustmap_input_options_block::PO_OUTPUT_TO_FILE
// clustmap_input_options_block::PO_SUMMARISE_TO_FILE

// eg_input_str();
// eg_input_mapfrom_str();
// eg_input_result_str();


BOOST_AUTO_TEST_CASE(processes_from_file_to_file) {
	// When calling perform_map_clusters with options: an input file, the output file option and an output file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_output_options_block::PO_OUTPUT_TO_FILE, TEMP_TEST_FILE_FILENAME.string() } );

	// Then expect:
	//  * an empty output stream and
	//  * the output file containing the correct output
	BOOST_CHECK_EQUAL( output_ss.str(), "" );
	BOOST_CHECK_FILES_EQUAL( TEMP_TEST_FILE_FILENAME, eg_renumber_only_result_file() );
}

// BOOST_AUTO_TEST_CASE(processes_from_stdin_to_output_file) {
// 	// Given an input stream containing the input
// 	input_ss.str( example_input_raw );

// 	// When calling perform_map_clusters with options: - (a dash to read from the input stream), the output file option and an output file
// 	execute_perform_map_clusters( {
// 		"-",
// 		"--" + clustmap_output_options_block::PO_QUIET,
// 		"--" + clustmap_output_options_block::PO_HITS_TEXT_TO_FILE, TEMP_TEST_FILE_FILENAME.string() } );

// 	// Then expect:
// 	//  * an empty output stream and
// 	//  * the output file containing the correct output
// 	//
// 	// \todo Add a better test tool for comparing a got file to an expected string
// 	BOOST_CHECK_EQUAL( output_ss.str(), "" );
// 	istringstream expected_out_ss{ example_output };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( expected_out_ss, "expected_crh_output", TEMP_TEST_FILE_FILENAME );
// }

// BOOST_AUTO_TEST_CASE(processes_from_stdin_to_output_file__deprecated_opts) {
// 	// Redirect any logging to log_ss
// 	log_to_ostream_guard log_output_guard{ log_ss };

// 	// Given an input stream containing the input
// 	input_ss.str( example_input_raw );

// 	// When calling perform_map_clusters with options: - (a dash to read from the input stream), the output file option and an output file
// 	execute_perform_map_clusters( { "-", "--" + clust_mapping_options_block::PO_OUTPUT_FILE, TEMP_TEST_FILE_FILENAME.string() } );

// 	// Then expect:
// 	//  * an empty output stream and
// 	//  * the output file containing the correct output
// 	//
// 	// \todo Add a better test tool for comparing a got file to an expected string
// 	BOOST_CHECK_EQUAL( output_ss.str(), "" );
// 	istringstream expected_out_ss{ example_output };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( expected_out_ss, "expected_crh_output", TEMP_TEST_FILE_FILENAME );

// 	BOOST_CHECK( regex_search( log_ss.str(), regex{ R"(deprecated.* \-\-quiet \-\-hits\-text\-to\-file )" } ) );
// }

// BOOST_AUTO_TEST_CASE(does_not_require_right_intersperses_all_to_cache) {
// 	// Given an input that requires caching at the start of match_b whilst processing in match_a and match_c
// 	// even though match_b doesn't right intersperse match_c
// 	const string input_hits_str =
// 		"query match_c 1 0-9,60-69\n"
// 		"query match_a 1 10-19,40-49\n"
// 		"query match_b 1 30-39,50-59\n";
// 	input_ss.str( input_hits_str );

// 	// When calling perform_map_clusters on that data with no trimming
// 	execute_perform_map_clusters( {
// 		"-",
// 		"--" + clust_mapping_options_block::PO_OVERLAP_TRIM_SPEC,
// 		"1/0"
// 	} );

// 	// Then expect the output to be the same as the input
// 	// (but with repeat of the boundaries for the resolved version)
// 	const string output_hits_str =
// 		"# Generated by cath-resolve-hits, one of the cath-tools (https://github.com/UCLOrengoGroup/cath-tools)\n"
// 		"#FIELDS query-id match-id score boundaries resolved\n"
// 		"query match_c 1 0-9,60-69 0-9,60-69\n"
// 		"query match_a 1 10-19,40-49 10-19,40-49\n"
// 		"query match_b 1 30-39,50-59 30-39,50-59\n";
// 	BOOST_CHECK_EQUAL( output_ss.str(), output_hits_str );
// }

// BOOST_AUTO_TEST_CASE(file_domtbl) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_DOMTBL_IN_FILENAME().string(), "--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_DOMTBL_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_DOMTBL_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_hmmsearch) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_hmmsearch_big_gap) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clustmap_input_options_block::PO_MIN_GAP_LENGTH, "10000"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_BIG_GAP_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_BIG_GAP_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_hmmsearch_small_gap) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clustmap_input_options_block::PO_MIN_GAP_LENGTH, "1"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_SMALL_GAP_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_SMALL_GAP_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_hmmsearch_trimmed) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clustmap_output_options_block::PO_OUTPUT_TRIMMED_HITS
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_TRIMMED_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_TRIMMED_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_hmmsearch_big_trim) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clust_mapping_options_block::PO_OVERLAP_TRIM_SPEC, "100/60"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_BIG_TRIM_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_BIG_TRIM_OUT_FILENAME() );
// }


// BOOST_AUTO_TEST_CASE(handles_output_hmmsearch_aln) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clustmap_output_options_block::PO_OUTPUT_HMMSEARCH_ALN
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_HMMSEARCH_ALN_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_HMMSEARCH_ALN_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_raw_evalue) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_RAW_EVALUE_IN_FILENAME().string(), "--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_EVALUES ),
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_RAW_EVALUE_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_RAW_EVALUE_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_raw_score) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_RAW_SCORE_IN_FILENAME().string(), "--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_SCORES ),
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_RAW_SCORE_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_RAW_SCORE_OUT_FILENAME() );
// }


// BOOST_AUTO_TEST_CASE(handles_dc_correctly) {
// 	execute_perform_map_clusters( {
// 		(CRH_CATH_DC_HANDLING_DATA_DIR() / "dc_eg_domtblout.in" ).string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
// 		"--" + crh_score_options_block::PO_APPLY_CATH_RULES,
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_CATH_DC_HANDLING_DATA_DIR() / "dc_eg_domtblout.cath_rules.out" );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_CATH_DC_HANDLING_DATA_DIR() / "dc_eg_domtblout.cath_rules.out" );
// }

// BOOST_AUTO_TEST_CASE(rejects_output_hmmsearch_aln_for_non_hmmsearch_format) {
// 	execute_perform_map_clusters( {
// 		(CRH_CATH_DC_HANDLING_DATA_DIR() / "dc_eg_domtblout.in" ).string(),
// 		"--" + clustmap_output_options_block::PO_OUTPUT_HMMSEARCH_ALN,
// 	} );
// 	BOOST_CHECK( boost::algorithm::contains( output_ss.str(), "Cannot use" ) );
// }

// BOOST_AUTO_TEST_CASE(generates_html_even_if_hmmsearch_aln_data_has_negative_scores) {
// 	const log_to_ostream_guard the_guard{ log_ss };

// 	BOOST_CHECK_NO_THROW_DIAG( execute_perform_map_clusters( {
// 		(CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.negatives_scores.in" ).string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clustmap_output_options_block::PO_HTML_OUTPUT_TO_FILE, "-"
// 	} ) );

// 	BOOST_CHECK( regex_search( log_ss.str(), regex{ R"(^Skipping .* weak hits\.)" } ) );
// }

// BOOST_AUTO_TEST_CASE(generates_html_even_if_hmmsearch_aln_data_has_negative_scores__deprecated_opts) {
// 	// Redirect any logging to log_ss
// 	log_to_ostream_guard log_output_guard{ log_ss };

// 	BOOST_CHECK_NO_THROW_DIAG( execute_perform_map_clusters( {
// 		(CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.negatives_scores.in" ).string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clust_mapping_options_block::PO_GENERATE_HTML_OUTPUT
// 	} ) );

// 	BOOST_CHECK( regex_search( log_ss.str(), regex{ R"(deprecated.* \-\-html\-output\-to\-file \-)" } ) );
// }

// BOOST_AUTO_TEST_CASE(handles_overlap_being_valid_only_once_short_seg_gone) {
// 	// Given: an input stream containing two hits that can both be included in the results
// 	// but only because the second hit's second segment is removed due to being shorter than
// 	// the min_seg_length (and not because of trimming)
// 	const string input_hits_str =
// 		"the_query match_1 1 21-122\n"
// 		"the_query match_2 1 0-20,94-95\n";
// 	input_ss.str( input_hits_str );

// 	// When: calling perform_map_clusters with options to read from the input stream,
// 	// remove segments shorter than 7 and perform no trimming
// 	execute_perform_map_clusters( {
// 		"-",
// 		"--" + clust_mapping_options_block::PO_MIN_SEG_LENGTH,    "7",
// 		"--" + clust_mapping_options_block::PO_OVERLAP_TRIM_SPEC, "1/0",
// 	} );

// 	// Then: expect the results to have included both hits but to have removed
// 	// the short segment from the final results
// 	const string output_hits_str =
// 		"# Generated by cath-resolve-hits, one of the cath-tools (https://github.com/UCLOrengoGroup/cath-tools)\n"
// 		"#FIELDS query-id match-id score boundaries resolved\n"
// 		"the_query match_2 1 0-20,94-95 0-20\n"
// 		"the_query match_1 1 21-122 21-122\n";
// 	BOOST_CHECK_EQUAL( output_ss.str(), output_hits_str );
// }

// BOOST_AUTO_TEST_SUITE(limit)

// BOOST_AUTO_TEST_CASE(file_domtbl) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_DOMTBL_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
// 		"--" + clust_mapping_options_block::PO_LIMIT_QUERIES + "=2"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_DOMTBL_LIMIT_2_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_DOMTBL_LIMIT_2_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_hmmsearch) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clust_mapping_options_block::PO_LIMIT_QUERIES + "=2"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_LIMIT_2_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_LIMIT_2_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_raw_evalue) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_RAW_EVALUE_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_EVALUES ),
// 		"--" + clust_mapping_options_block::PO_LIMIT_QUERIES + "=2"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_RAW_EVALUE_LIMIT_2_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_RAW_EVALUE_LIMIT_2_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(file_raw_score) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_RAW_SCORE_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_SCORES ),
// 		"--" + clust_mapping_options_block::PO_LIMIT_QUERIES + "=2"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_RAW_SCORE_LIMIT_2_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_RAW_SCORE_LIMIT_2_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_SUITE_END()



// BOOST_AUTO_TEST_SUITE(summary_output)

// BOOST_AUTO_TEST_CASE(summarise) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clustmap_output_options_block::PO_SUMMARISE_TO_FILE, "-"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_SUMMARISE_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_SUMMARISE_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(summarise__deprecated_opts) {
// 	// Redirect any logging to log_ss
// 	log_to_ostream_guard log_output_guard{ log_ss };

// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clust_mapping_options_block::PO_SUMMARISE
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_SUMMARISE_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_SUMMARISE_OUT_FILENAME() );

// 	BOOST_CHECK( regex_search( log_ss.str(), regex{ R"(deprecated.* \-\-summarise\-to\-file \-)" } ) );
// }

// BOOST_AUTO_TEST_SUITE_END()


// BOOST_AUTO_TEST_SUITE(css)

// BOOST_AUTO_TEST_CASE(export_css_to_stdout) {
// 	execute_perform_map_clusters( {
// 		"--" + clustmap_output_options_block::PO_EXPORT_CSS_FILE, "-",
// 	} );

// 	BOOST_CHECK(   regex_search( output_ss.str(), regex{ R"(^/\* \-\-\- Start[^]*crh-exclusion-note[^]{0,100}$)" } ) ); // In EMCAScript, . doesn't match newline characters, hence the use of [^] here
// 	BOOST_CHECK( ! regex_search( output_ss.str(), regex{ R"(443cb81e8e280e529de69ef113974208)" } ) );
// }

// BOOST_AUTO_TEST_CASE(export_css_to_file) {
// 	execute_perform_map_clusters( {
// 		"--" + clustmap_output_options_block::PO_EXPORT_CSS_FILE, TEMP_TEST_FILE_FILENAME.string(),
// 	} );

// 	BOOST_CHECK(   output_ss.str().empty() );
// 	BOOST_CHECK(   regex_search( read_string_from_file( TEMP_TEST_FILE_FILENAME ), regex{ R"(^/\* \-\-\- Start[^]*crh-exclusion-note[^]{0,100}$)" } ) ); // In EMCAScript, . doesn't match newline characters, hence the use of [^] here
// 	BOOST_CHECK( ! regex_search( read_string_from_file( TEMP_TEST_FILE_FILENAME ), regex{ R"(443cb81e8e280e529de69ef113974208)" } ) );
// }

// BOOST_AUTO_TEST_SUITE_END()


// BOOST_AUTO_TEST_SUITE(html_output)

// BOOST_AUTO_TEST_CASE(html) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clust_mapping_options_block::PO_OVERLAP_TRIM_SPEC, "150/90",
// 		"--" + clustmap_output_options_block::PO_HTML_OUTPUT_TO_FILE, "-",
// 		"--" + clust_mapping_options_block::PO_WORST_PERMISSIBLE_BITSCORE, "14",
// 		"--" + clust_mapping_options_block::PO_EXCLUDE_REJECTED_HITS,
// 		"--" + clust_mapping_options_block::PO_MAX_NUM_NON_SOLN_HITS, "10"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_HTML_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_HTML_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(html__deprecated_opts) {
// 	// Redirect any logging to log_ss
// 	log_to_ostream_guard log_output_guard{ log_ss };

// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clust_mapping_options_block::PO_OVERLAP_TRIM_SPEC, "150/90",
// 		"--" + clust_mapping_options_block::PO_GENERATE_HTML_OUTPUT,
// 		"--" + clust_mapping_options_block::PO_WORST_PERMISSIBLE_BITSCORE, "14",
// 		"--" + clust_mapping_options_block::PO_EXCLUDE_REJECTED_HITS,
// 		"--" + clust_mapping_options_block::PO_MAX_NUM_NON_SOLN_HITS, "10"
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_HTML_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_HTML_OUT_FILENAME() );

// 	BOOST_CHECK( regex_search( log_ss.str(), regex{ R"(deprecated.* \-\-html\-output\-to\-file \-)" } ) );
// }

// BOOST_AUTO_TEST_SUITE_END()


// BOOST_AUTO_TEST_SUITE(json_output)

// BOOST_AUTO_TEST_CASE(json_from_hmmsearch_out) {
// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clustmap_output_options_block::PO_JSON_OUTPUT_TO_FILE, "-",
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_JSON_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_JSON_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(json_from_domtblout) {
// 	execute_perform_map_clusters( {
// 		(CRH_TEST_DATA_DIR() / "eg_domtblout.in" ).string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
// 		"--" + clustmap_output_options_block::PO_JSON_OUTPUT_TO_FILE, "-",
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_DOMTBL_JSON_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_DOMTBL_JSON_OUT_FILENAME() );
// }

// BOOST_AUTO_TEST_CASE(json_from_hmmsearch_out__deprecated_opts) {
// 	// Redirect any logging to log_ss
// 	log_to_ostream_guard log_output_guard{ log_ss };

// 	execute_perform_map_clusters( {
// 		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
// 		"--" + clust_mapping_options_block::PO_JSON_OUTPUT,
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_JSON_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_JSON_OUT_FILENAME() );

// 	BOOST_CHECK( regex_search( log_ss.str(), regex{ R"(deprecated.* \-\-json\-output\-to\-file \-)" } ) );
// }

// BOOST_AUTO_TEST_CASE(json_from_domtblout__deprecated_opts) {
// 	// Redirect any logging to log_ss
// 	log_to_ostream_guard log_output_guard{ log_ss };

// 	execute_perform_map_clusters( {
// 		(CRH_TEST_DATA_DIR() / "eg_domtblout.in" ).string(),
// 		"--" + clustmap_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
// 		"--" + clust_mapping_options_block::PO_JSON_OUTPUT,
// 	} );
// 	istringstream istream_of_output{ output_ss.str() };
// 	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_DOMTBL_JSON_OUT_FILENAME() );
// 	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_DOMTBL_JSON_OUT_FILENAME() );

// 	BOOST_CHECK( regex_search( log_ss.str(), regex{ R"(deprecated.* \-\-json\-output\-to\-file \-)" } ) );
// }

// BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()
