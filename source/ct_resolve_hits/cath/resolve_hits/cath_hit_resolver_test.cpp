/// \file
/// \brief The cath_hit_resolver test suite

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

#include <filesystem>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/join.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/boost_addenda/log/stringstream_log_sink.hpp"
#include "cath/common/file/read_string_from_file.hpp"
#include "cath/common/file/simple_file_read_write.hpp"
#include "cath/common/file/temp_file.hpp"
#include "cath/resolve_hits/cath_hit_resolver.hpp"
#include "cath/resolve_hits/options/options_block/crh_filter_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_html_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_input_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_output_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_score_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_segment_options_block.hpp"
#include "cath/resolve_hits/options/options_block/crh_single_output_options_block.hpp"
#include "cath/resolve_hits/test/resolve_hits_fixture.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "cath/test/global_test_constants.hpp"
#include "cath/test/predicate/files_equal.hpp"
#include "cath/test/predicate/string_matches_file.hpp"

#include <regex>

namespace cath { namespace test { } }

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::cath::rslv;
using namespace ::cath::test;
using namespace ::std::literals::string_literals;

using ::boost::algorithm::contains;
using ::boost::range::join;
using ::cath::common::copy_build;
using ::cath::common::temp_file;
using ::cath::common::write_file;
using ::std::filesystem::path;
using ::std::istringstream;
using ::std::ostringstream;
using ::std::regex;
using ::std::string;

namespace cath {
	namespace test {

		/// \brief The hit_resolver_test_suite_fixture to assist in testing calc_hit_list
		struct hit_resolver_test_suite_fixture : protected resolve_hits_fixture, protected global_test_constants {
		protected:
			~hit_resolver_test_suite_fixture() noexcept  = default;

			/// \brief Call perform_resolve_hits() with the specified arguments preceded by a pseudo-program-name
			///        and the fixture's i/o streams.
			void execute_perform_resolve_hits(const str_vec &prm_arguments ///< The arguments to pass to perform_resolve_hits(), preceded by a pseudo-program-name
			                                  ) {
				const auto progname_vec = { "pseudo_program_name"s };
				perform_resolve_hits(
					copy_build<str_vec>( join(
						progname_vec,
						prm_arguments
					) ),
					input_ss, output_ss, parse_sources::CMND_LINE_ONLY
				);
			}

			/// \brief The input stream to use in the tests
			istringstream   input_ss;
			
			/// \brief The output stream to use in the tests
			ostringstream   output_ss;

			/// \brief A temporary temp_file
			const temp_file TEMP_TEST_FILE{ ".cath_hit_resolver__test_file.%%%%-%%%%-%%%%-%%%%" };

			/// \brief The path of the TEMP_TEST_FILE
			const path      TEMP_TEST_FILE_FILENAME = get_filename( TEMP_TEST_FILE );
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(cath_hit_resolver_test_suite, hit_resolver_test_suite_fixture)

BOOST_AUTO_TEST_CASE(fails_on_attempt_to_mix_deprecated_options_with_new) {
	// When calling perform_resolve_hits with options: - (a dash to read from the input stream)
	execute_perform_resolve_hits( {
		"--" + crh_output_options_block::PO_JSON_OUTPUT_TO_FILE, "-"
		"--" + crh_single_output_options_block::PO_SUMMARISE,
		"--" + crh_single_output_options_block::PO_OUTPUT_FILE, TEMP_TEST_FILE_FILENAME.string(),
		"--" + crh_output_options_block::PO_QUIET,
	} );

	// Then expect the correct error message in the output stream
	BOOST_CHECK_EQUAL( output_ss.str(), "cath-resolve-hits: Cannot mix old, deprecated options (--output-file) with new, replacement options (--quiet, --json-output-to-file). Please use the new options only.\nSee 'cath-resolve-hits -h' for usage.\n" );
}


BOOST_AUTO_TEST_CASE(processes_from_stdin_to_stdout) {
	// Given an input stream containing the input
	input_ss.str( example_input_raw );

	// When calling perform_resolve_hits with options: - (a dash to read from the input stream)
	execute_perform_resolve_hits( { "-" } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_EQUAL( blank_vrsn( output_ss ), example_output );
}


BOOST_AUTO_TEST_CASE(processes_from_file_to_stdout) {
	// Given an input file containing the input data
	write_file( TEMP_TEST_FILE_FILENAME, example_input_raw );

	// When calling perform_resolve_hits with options: the input file
	execute_perform_resolve_hits( { TEMP_TEST_FILE_FILENAME.string() } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_EQUAL( blank_vrsn( output_ss ), example_output );
}

BOOST_AUTO_TEST_CASE(processes_from_stdin_to_output_file) {
	// Given an input stream containing the input
	input_ss.str( example_input_raw );

	// When calling perform_resolve_hits with options: - (a dash to read from the input stream), the output file option and an output file
	execute_perform_resolve_hits( {
		"-",
		"--" + crh_output_options_block::PO_QUIET,
		"--" + crh_output_options_block::PO_HITS_TEXT_TO_FILE, TEMP_TEST_FILE_FILENAME.string() } );

	// Then expect:
	//  * an empty output stream and
	//  * the output file containing the correct output
	//
	// \todo Add a better test tool for comparing a got file to an expected string
	BOOST_CHECK_EQUAL( output_ss.str(), "" );
	BOOST_CHECK_STRING_MATCHES_FILE( example_output, blank_vrsn( TEMP_TEST_FILE_FILENAME ) );
}

BOOST_AUTO_TEST_CASE(processes_from_stdin_to_output_file__deprecated_opts) {
	// Redirect any logging
	stringstream_log_sink log_sink;

	// Given an input stream containing the input
	input_ss.str( example_input_raw );

	// When calling perform_resolve_hits with options: - (a dash to read from the input stream), the output file option and an output file
	execute_perform_resolve_hits( { "-", "--" + crh_single_output_options_block::PO_OUTPUT_FILE, TEMP_TEST_FILE_FILENAME.string() } );

	// Then expect:
	//  * an empty output stream and
	//  * the output file containing the correct output
	//
	// \todo Add a better test tool for comparing a got file to an expected string
	BOOST_CHECK_EQUAL( output_ss.str(), "" );
	BOOST_CHECK_STRING_MATCHES_FILE( example_output, blank_vrsn( TEMP_TEST_FILE_FILENAME ) );

	BOOST_CHECK( regex_search( log_sink.str(), regex{ R"(deprecated.* \-\-quiet \-\-hits\-text\-to\-file )" } ) );
}

BOOST_AUTO_TEST_CASE(does_not_require_right_intersperses_all_to_cache) {
	// Given an input that requires caching at the start of match_b whilst processing in match_a and match_c
	// even though match_b doesn't right intersperse match_c
	const string input_hits_str =
		"query match_c 1 0-9,60-69\n"
		"query match_a 1 10-19,40-49\n"
		"query match_b 1 30-39,50-59\n";
	input_ss.str( input_hits_str );

	// When calling perform_resolve_hits on that data with no trimming
	execute_perform_resolve_hits( {
		"-",
		"--" + crh_segment_options_block::PO_OVERLAP_TRIM_SPEC,
		"1/0"
	} );
	
	// Then expect the output to be the same as the input
	// (but with repeat of the boundaries for the resolved version)
	const string output_hits_str =
		"# Generated by cath-resolve-hits vX.X.X-X-XXXXXXXX, one of the cath-tools (https://github.com/UCLOrengoGroup/cath-tools)\n"
		"#FIELDS query-id match-id score boundaries resolved\n"
		"query match_c 1 0-9,60-69 0-9,60-69\n"
		"query match_a 1 10-19,40-49 10-19,40-49\n"
		"query match_b 1 30-39,50-59 30-39,50-59\n";
	BOOST_CHECK_EQUAL( blank_vrsn( output_ss ), output_hits_str );
}

BOOST_AUTO_TEST_CASE(file_domtbl) {
	execute_perform_resolve_hits( {
		CRH_EG_DOMTBL_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_DOMTBL_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_HMMSEARCH_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch_big_gap) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_input_options_block::PO_MIN_GAP_LENGTH, "10000"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_HMMSEARCH_BIG_GAP_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch_small_gap) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_input_options_block::PO_MIN_GAP_LENGTH, "1"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_HMMSEARCH_SMALL_GAP_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch_trimmed) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_OUTPUT_TRIMMED_HITS
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_HMMSEARCH_TRIMMED_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch_big_trim) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_segment_options_block::PO_OVERLAP_TRIM_SPEC, "100/60"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_HMMSEARCH_BIG_TRIM_OUT_FILENAME() );
}


BOOST_AUTO_TEST_CASE(handles_output_hmmer_aln) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_OUTPUT_HMMER_ALN
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_HMMSEARCH_HMMSEARCH_ALN_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_raw_evalue) {
	execute_perform_resolve_hits( {
		CRH_EG_RAW_EVALUE_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_EVALUES ),
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_RAW_EVALUE_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_raw_score) {
	execute_perform_resolve_hits( {
		CRH_EG_RAW_SCORE_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_SCORES ),
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_RAW_SCORE_OUT_FILENAME() );
}


BOOST_AUTO_TEST_CASE(handles_dc_correctly) {
	execute_perform_resolve_hits( {
		(CRH_CATH_DC_HANDLING_DATA_DIR() / "dc_eg_domtblout.in" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
		"--" + crh_score_options_block::PO_APPLY_CATH_RULES,
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_CATH_DC_HANDLING_DATA_DIR() / "dc_eg_domtblout.cath_rules.out" );
}

BOOST_AUTO_TEST_CASE(rejects_output_hmmer_aln_for_non_hmmsearch_format) {
	execute_perform_resolve_hits( {
		(CRH_CATH_DC_HANDLING_DATA_DIR() / "dc_eg_domtblout.in" ).string(),
		"--" + crh_output_options_block::PO_OUTPUT_HMMER_ALN,
	} );
	BOOST_CHECK( boost::algorithm::contains( output_ss.str(), "Cannot use" ) );
}

BOOST_AUTO_TEST_CASE(generates_html_even_if_hmmsearch_aln_data_has_negative_scores) {
	const stringstream_log_sink log_sink;

	BOOST_CHECK_NO_THROW_DIAG( execute_perform_resolve_hits( {
		(CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.negatives_scores.in" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_HTML_OUTPUT_TO_FILE, "-"
	} ) );

	BOOST_CHECK( regex_search( log_sink.str(), regex{ R"(^Skipping .* weak hits\.)" } ) );
}

BOOST_AUTO_TEST_CASE(generates_html_even_if_hmmsearch_aln_data_has_negative_scores__deprecated_opts) {
	// Redirect any logging
	stringstream_log_sink log_sink;

	BOOST_CHECK_NO_THROW_DIAG( execute_perform_resolve_hits( {
		(CRH_TEST_DATA_DIR() / "eg_hmmsearch_out.negatives_scores.in" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_single_output_options_block::PO_GENERATE_HTML_OUTPUT
	} ) );

	BOOST_CHECK( regex_search( log_sink.str(), regex{ R"(deprecated.* \-\-html\-output\-to\-file \-)" } ) );
}

BOOST_AUTO_TEST_CASE(handles_overlap_being_valid_only_once_short_seg_gone) {
	// Given: an input stream containing two hits that can both be included in the results
	// but only because the second hit's second segment is removed due to being shorter than
	// the min_seg_length (and not because of trimming)
	const string input_hits_str =
		"the_query match_1 1 21-122\n"
		"the_query match_2 1 0-20,94-95\n";
	input_ss.str( input_hits_str );

	// When: calling perform_resolve_hits with options to read from the input stream,
	// remove segments shorter than 7 and perform no trimming
	execute_perform_resolve_hits( {
		"-",
		"--" + crh_segment_options_block::PO_MIN_SEG_LENGTH,    "7",
		"--" + crh_segment_options_block::PO_OVERLAP_TRIM_SPEC, "1/0",
	} );

	// Then: expect the results to have included both hits but to have removed
	// the short segment from the final results
	const string output_hits_str =
		"# Generated by cath-resolve-hits vX.X.X-X-XXXXXXXX, one of the cath-tools (https://github.com/UCLOrengoGroup/cath-tools)\n"
		"#FIELDS query-id match-id score boundaries resolved\n"
		"the_query match_2 1 0-20,94-95 0-20\n"
		"the_query match_1 1 21-122 21-122\n";
	BOOST_CHECK_EQUAL( blank_vrsn( output_ss ), output_hits_str );
}

BOOST_AUTO_TEST_SUITE(limit)

BOOST_AUTO_TEST_CASE(file_domtbl) {
	execute_perform_resolve_hits( {
		CRH_EG_DOMTBL_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
		"--" + crh_filter_options_block::PO_LIMIT_QUERIES + "=2"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_DOMTBL_LIMIT_2_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_filter_options_block::PO_LIMIT_QUERIES + "=2"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_HMMSEARCH_LIMIT_2_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_raw_evalue) {
	execute_perform_resolve_hits( {
		CRH_EG_RAW_EVALUE_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_EVALUES ),
		"--" + crh_filter_options_block::PO_LIMIT_QUERIES + "=2"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_RAW_EVALUE_LIMIT_2_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_raw_score) {
	execute_perform_resolve_hits( {
		CRH_EG_RAW_SCORE_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_SCORES ),
		"--" + crh_filter_options_block::PO_LIMIT_QUERIES + "=2"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( blank_vrsn( output_ss ), CRH_EG_RAW_SCORE_LIMIT_2_OUT_FILENAME() );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(summary_output)

BOOST_AUTO_TEST_CASE(summarise) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_SUMMARISE_TO_FILE, "-"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), CRH_EG_HMMSEARCH_SUMMARISE_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(summarise__deprecated_opts) {
	// Redirect any logging
	stringstream_log_sink log_sink;

	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_single_output_options_block::PO_SUMMARISE
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), CRH_EG_HMMSEARCH_SUMMARISE_OUT_FILENAME() );

	BOOST_CHECK( regex_search( log_sink.str(), regex{ R"(deprecated.* \-\-summarise\-to\-file \-)" } ) );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(css)

BOOST_AUTO_TEST_CASE(export_css_to_stdout) {
	execute_perform_resolve_hits( {
		"--" + crh_output_options_block::PO_EXPORT_CSS_FILE, "-",
	} );

	BOOST_CHECK(   regex_search( output_ss.str(), regex{ R"(^/\* \-\-\- Start[^]*crh-exclusion-note[^]{0,100}$)" } ) ); // In EMCAScript, . doesn't match newline characters, hence the use of [^] here
	BOOST_CHECK( ! regex_search( output_ss.str(), regex{ R"(443cb81e8e280e529de69ef113974208)" } ) );
}

BOOST_AUTO_TEST_CASE(export_css_to_file) {
	execute_perform_resolve_hits( {
		"--" + crh_output_options_block::PO_EXPORT_CSS_FILE, TEMP_TEST_FILE_FILENAME.string(),
	} );

	BOOST_CHECK(   output_ss.str().empty() );
	BOOST_CHECK(   regex_search( read_string_from_file( TEMP_TEST_FILE_FILENAME ), regex{ R"(^/\* \-\-\- Start[^]*crh-exclusion-note[^]{0,100}$)" } ) ); // In EMCAScript, . doesn't match newline characters, hence the use of [^] here
	BOOST_CHECK( ! regex_search( read_string_from_file( TEMP_TEST_FILE_FILENAME ), regex{ R"(443cb81e8e280e529de69ef113974208)" } ) );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(html_output)

BOOST_AUTO_TEST_CASE(html) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_segment_options_block::PO_OVERLAP_TRIM_SPEC, "150/90",
		"--" + crh_output_options_block::PO_HTML_OUTPUT_TO_FILE, "-",
		"--" + crh_filter_options_block::PO_WORST_PERMISSIBLE_BITSCORE, "14",
		"--" + crh_html_options_block::PO_EXCLUDE_REJECTED_HITS,
		"--" + crh_html_options_block::PO_MAX_NUM_NON_SOLN_HITS, "10"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), CRH_EG_HMMSEARCH_HTML_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(html__deprecated_opts) {
	// Redirect any logging
	stringstream_log_sink log_sink;

	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_segment_options_block::PO_OVERLAP_TRIM_SPEC, "150/90",
		"--" + crh_single_output_options_block::PO_GENERATE_HTML_OUTPUT,
		"--" + crh_filter_options_block::PO_WORST_PERMISSIBLE_BITSCORE, "14",
		"--" + crh_html_options_block::PO_EXCLUDE_REJECTED_HITS,
		"--" + crh_html_options_block::PO_MAX_NUM_NON_SOLN_HITS, "10"
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), CRH_EG_HMMSEARCH_HTML_OUT_FILENAME() );

	BOOST_CHECK( regex_search( log_sink.str(), regex{ R"(deprecated.* \-\-html\-output\-to\-file \-)" } ) );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(json_output)

BOOST_AUTO_TEST_CASE(json_from_hmmsearch_out) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_JSON_OUTPUT_TO_FILE, "-",
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), CRH_EG_HMMSEARCH_JSON_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(json_from_domtblout) {
	execute_perform_resolve_hits( {
		(CRH_TEST_DATA_DIR() / "eg_domtblout.in" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
		"--" + crh_output_options_block::PO_JSON_OUTPUT_TO_FILE, "-",
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), CRH_EG_DOMTBL_JSON_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(json_from_hmmsearch_out__deprecated_opts) {
	// Redirect any logging
	stringstream_log_sink log_sink;

	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_single_output_options_block::PO_JSON_OUTPUT,
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), CRH_EG_HMMSEARCH_JSON_OUT_FILENAME() );

	BOOST_CHECK( regex_search( log_sink.str(), regex{ R"(deprecated.* \-\-json\-output\-to\-file \-)" } ) );
}

BOOST_AUTO_TEST_CASE(json_from_domtblout__deprecated_opts) {
	// Redirect any logging
	stringstream_log_sink log_sink;

	execute_perform_resolve_hits( {
		(CRH_TEST_DATA_DIR() / "eg_domtblout.in" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTBLOUT ),
		"--" + crh_single_output_options_block::PO_JSON_OUTPUT,
	} );
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), CRH_EG_DOMTBL_JSON_OUT_FILENAME() );

	BOOST_CHECK( regex_search( log_sink.str(), regex{ R"(deprecated.* \-\-json\-output\-to\-file \-)" } ) );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(hmm_coverage)

BOOST_AUTO_TEST_CASE(hmm_coverage__neither) {
	execute_perform_resolve_hits( {
		(CRH_HMM_COVERAGE_DATA_DIR() / "hmm_coverage.hmmsearch_out.in" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT,         to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_HITS_TEXT_TO_FILE,   TEMP_TEST_FILE_FILENAME.string()
	} );
	BOOST_CHECK_FILES_EQUAL( blank_vrsn( TEMP_TEST_FILE_FILENAME ), CRH_HMM_COVERAGE_DATA_DIR() / "hmm_coverage.out.cov_neither" );
}

BOOST_AUTO_TEST_CASE(hmm_coverage__discontinuous) {
	execute_perform_resolve_hits( {
		(CRH_HMM_COVERAGE_DATA_DIR() / "hmm_coverage.hmmsearch_out.in" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT,       to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_HITS_TEXT_TO_FILE, TEMP_TEST_FILE_FILENAME.string(),
		"--" + crh_filter_options_block::PO_MIN_DC_HMM_COVERAGE + "=99.9"
	} );
	BOOST_CHECK_FILES_EQUAL( blank_vrsn( TEMP_TEST_FILE_FILENAME ), CRH_HMM_COVERAGE_DATA_DIR() / "hmm_coverage.out.cov_dc" );
}

BOOST_AUTO_TEST_CASE(hmm_coverage__all) {
	execute_perform_resolve_hits( {
		(CRH_HMM_COVERAGE_DATA_DIR() / "hmm_coverage.hmmsearch_out.in" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT,       to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_HITS_TEXT_TO_FILE, TEMP_TEST_FILE_FILENAME.string(),
		"--" + crh_filter_options_block::PO_MIN_HMM_COVERAGE    + "=60.0"
	} );
	BOOST_CHECK_FILES_EQUAL( blank_vrsn( TEMP_TEST_FILE_FILENAME ), CRH_HMM_COVERAGE_DATA_DIR() / "hmm_coverage.out.cov_all" );
}

BOOST_AUTO_TEST_CASE(hmm_coverage__all_and_discontinuous) {
	execute_perform_resolve_hits( {
		(CRH_HMM_COVERAGE_DATA_DIR() / "hmm_coverage.hmmsearch_out.in" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT,       to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_HITS_TEXT_TO_FILE, TEMP_TEST_FILE_FILENAME.string(),
		"--" + crh_filter_options_block::PO_MIN_DC_HMM_COVERAGE + "=99.9",
		"--" + crh_filter_options_block::PO_MIN_HMM_COVERAGE    + "=60.0"
	} );
	BOOST_CHECK_FILES_EQUAL( blank_vrsn( TEMP_TEST_FILE_FILENAME ), CRH_HMM_COVERAGE_DATA_DIR() / "hmm_coverage.out.cov_both" );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(hmmscan_format)

BOOST_AUTO_TEST_CASE(seqs_hmmsearch) {
	execute_perform_resolve_hits( {
		( CRH_HMMSCAN_DATA_DIR() / "seqs.hmmsearch" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT,       to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_HITS_TEXT_TO_FILE, TEMP_TEST_FILE_FILENAME.string()
	} );
	BOOST_CHECK_FILES_EQUAL( blank_vrsn( TEMP_TEST_FILE_FILENAME ), CRH_HMMSCAN_DATA_DIR() / "seqs.hmmsearch.out" );
}

BOOST_AUTO_TEST_CASE(seqs_hmmscan) {
	execute_perform_resolve_hits( {
		( CRH_HMMSCAN_DATA_DIR() / "seqs.hmmscan" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT,       to_string( hits_input_format_tag::HMMSCAN_OUT ),
		"--" + crh_output_options_block::PO_HITS_TEXT_TO_FILE, TEMP_TEST_FILE_FILENAME.string()
	} );
	BOOST_CHECK_FILES_EQUAL( blank_vrsn( TEMP_TEST_FILE_FILENAME ), CRH_HMMSCAN_DATA_DIR() / "seqs.hmmscan.out" );
}

BOOST_AUTO_TEST_CASE(p53_p63_hmmscan) {
	execute_perform_resolve_hits( {
		( CRH_HMMSCAN_DATA_DIR() / "p53_p63.hmmscan" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT,       to_string( hits_input_format_tag::HMMSCAN_OUT ),
		"--" + crh_output_options_block::PO_HITS_TEXT_TO_FILE, TEMP_TEST_FILE_FILENAME.string()
	} );
	BOOST_CHECK_FILES_EQUAL( blank_vrsn( TEMP_TEST_FILE_FILENAME ), CRH_HMMSCAN_DATA_DIR() / "p53_p63.hmmscan.out" );
}

BOOST_AUTO_TEST_CASE(single_sequence_hmmscan_out) {
	execute_perform_resolve_hits( {
		( CRH_HMMSCAN_DATA_DIR() / "single_sequence.hmmscan.out" ).string(),
		"--" + crh_input_options_block::PO_INPUT_FORMAT,       to_string( hits_input_format_tag::HMMSCAN_OUT ),
		"--" + crh_output_options_block::PO_HITS_TEXT_TO_FILE, TEMP_TEST_FILE_FILENAME.string()
	} );
	BOOST_CHECK_FILES_EQUAL( blank_vrsn( TEMP_TEST_FILE_FILENAME ), CRH_HMMSCAN_DATA_DIR() / "single_sequence.hmmscan.out.out" );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
