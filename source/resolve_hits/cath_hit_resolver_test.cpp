
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

#include <boost/range/join.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/algorithm/copy_build.h"
#include "common/file/simple_file_read_write.h"
#include "common/file/temp_file.h"
#include "common/test_predicate/istream_and_file_equal.h"
#include "resolve_hits/cath_hit_resolver.h"
#include "resolve_hits/options/options_block/crh_input_options_block.h"
#include "resolve_hits/options/options_block/crh_output_options_block.h"
#include "resolve_hits/options/options_block/crh_segment_options_block.h"
#include "test/global_test_constants.h"
#include "test/resolve_hits/resolve_hits_fixture.h"

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace cath::test;

using boost::filesystem::path;
using boost::range::join;
using std::istringstream;
using std::ostringstream;
using std::string;

namespace cath {
	namespace test {

		/// \brief The hit_resolver_test_suite_fixture to assist in testing calc_hit_list
		struct hit_resolver_test_suite_fixture : protected resolve_hits_fixture, protected global_test_constants {
		protected:
			~hit_resolver_test_suite_fixture() noexcept  = default;

			/// \brief Call perform_resolve_hits() with the specified arguments preceded by a pseudo-program-name
			///        and the fixture's i/o streams. Test that the error stream is empty afterwards.
			void execute_perform_resolve_hits(const str_vec &arg_arguments ///< The arguments to pass to perform_resolve_hits(), preceded by a pseudo-program-name
			                                  ) {
				const str_vec progname_vec = { "pseudo_program_name" };
				perform_resolve_hits(
					copy_build<str_vec>( join(
						progname_vec,
						arg_arguments
					) ),
					input_ss, output_ss, error_ss
				);

				// Tests for crh options errors should be in chr_options_test.cpp, not here
				// so ensure that the errror error_ss is empty at the end of every perform_resolve_hits()
				BOOST_CHECK_EQUAL( error_ss.str(), "" );
			}

			/// \brief The input stream to use in the tests
			istringstream   input_ss;
			
			/// \brief The output stream to use in the tests
			ostringstream   output_ss;

			/// \brief The error stream to use in the tests
			ostringstream   error_ss;

			/// \brief A temporary temp_file
			const temp_file TEMP_TEST_FILE{ ".cath_hit_resolver__test_file.%%%%-%%%%-%%%%-%%%%" };

			/// \brief The path of the TEMP_TEST_FILE
			const path      TEMP_TEST_FILE_FILENAME = get_filename( TEMP_TEST_FILE );
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(cath_hit_resolver_test_suite, hit_resolver_test_suite_fixture)


BOOST_AUTO_TEST_CASE(processes_from_stdin_to_stdout) {
	// Given an input stream containing the input
	input_ss.str( example_input_raw );

	// When calling perform_resolve_hits with options: - (a dash to read from the input stream)
	execute_perform_resolve_hits( { "-" } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_EQUAL( output_ss.str(), example_output );
}


BOOST_AUTO_TEST_CASE(processes_from_file_to_stdout) {
	// Given an input file containing the input data
	write_file( TEMP_TEST_FILE_FILENAME, example_input_raw );

	// When calling perform_resolve_hits with options: the input file
	execute_perform_resolve_hits( { TEMP_TEST_FILE_FILENAME.string() } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_EQUAL( output_ss.str(), example_output );
}


BOOST_AUTO_TEST_CASE(processes_from_stdin_to_output_file) {
	// Given an input stream containing the input
	input_ss.str( example_input_raw );

	// When calling perform_resolve_hits with options: - (a dash to read from the input stream), the output file option and an output file
	execute_perform_resolve_hits( { "-", "--" + crh_output_options_block::PO_OUTPUT_FILE, TEMP_TEST_FILE_FILENAME.string() } );

	// Then expect:
	//  * an empty output stream and
	//  * the output file containing the correct output
	//
	// \todo Add a better test tool for comparing a got file to an expected string
	BOOST_CHECK_EQUAL( output_ss.str(), "" );
	istringstream expected_out_ss{ example_output };
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( expected_out_ss, "expected_crh_output", TEMP_TEST_FILE_FILENAME );
}

BOOST_AUTO_TEST_CASE(does_not_require_right_intersperses_all_to_cache) {
	// Given an input that requires caching at the start of match_b whilst processing in match_a and match_c
	// even though match_b doesn't right intersperse match_c
	const string hits_str =
		"query match_c 1 0-9,60-69\n"
		"query match_a 1 10-19,40-49\n"
		"query match_b 1 30-39,50-59\n";
	input_ss.str( hits_str );

	// When calling perform_resolve_hits on that data with no trimming
	execute_perform_resolve_hits( {
		"-",
		"--" + crh_segment_options_block::PO_OVERLAP_TRIM_SPEC,
		"1/0"
	} );
	
	// Then expect the output to be the same as the input
	BOOST_CHECK_EQUAL( output_ss.str(), hits_str );
}

BOOST_AUTO_TEST_CASE(file_domtbl) {
	execute_perform_resolve_hits( {
		CRH_EG_DOMTBL_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMER_DOMTMBLOUT ),
	} );
	istringstream istream_of_output{ output_ss.str() };
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_DOMTBL_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
	} );
	istringstream istream_of_output{ output_ss.str() };
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch_big_gap) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_input_options_block::PO_MIN_GAP_LENGTH, "10000"
	} );
	istringstream istream_of_output{ output_ss.str() };
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_BIG_GAP_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch_small_gap) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_input_options_block::PO_MIN_GAP_LENGTH, "1"
	} );
	istringstream istream_of_output{ output_ss.str() };
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_SMALL_GAP_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_hmmsearch_trimmed) {
	execute_perform_resolve_hits( {
		CRH_EG_HMMSEARCH_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::HMMSEARCH_OUT ),
		"--" + crh_output_options_block::PO_OUTPUT_TRIMMED_HITS
	} );
	istringstream istream_of_output{ output_ss.str() };
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_TRIMMED_OUT_FILENAME() );
	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_HMMSEARCH_TRIMMED_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_raw_evalue) {
	execute_perform_resolve_hits( {
		CRH_EG_RAW_EVALUE_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_EVALUES ),
	} );
	istringstream istream_of_output{ output_ss.str() };
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL             ( istream_of_output, "got_ss", CRH_EG_RAW_EVALUE_OUT_FILENAME() );
	// BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( istream_of_output, "got_ss", CRH_EG_RAW_EVALUE_OUT_FILENAME() );
}

BOOST_AUTO_TEST_CASE(file_raw_score) {
	execute_perform_resolve_hits( {
		CRH_EG_RAW_SCORE_IN_FILENAME().string(), "--" + crh_input_options_block::PO_INPUT_FORMAT, to_string( hits_input_format_tag::RAW_WITH_SCORES ),
	} );
	istringstream istream_of_output{ output_ss.str() };
	BOOST_CHECK_ISTREAM_AND_FILE_EQUAL( istream_of_output, "got_ss", CRH_EG_RAW_SCORE_OUT_FILENAME() );
}

BOOST_AUTO_TEST_SUITE_END()
