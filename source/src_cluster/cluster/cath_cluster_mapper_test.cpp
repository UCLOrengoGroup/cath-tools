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
#include <boost/test/unit_test.hpp>

#include "cluster/cath_cluster_mapper.hpp"
#include "cluster/options/options_block/clust_mapping_options_block.hpp"
#include "cluster/options/options_block/clustmap_input_options_block.hpp"
#include "cluster/options/options_block/clustmap_output_options_block.hpp"
#include "cluster/test/map_clusters_fixture.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/file/temp_file.hpp"
#include "common/regex/regex_count.hpp"
#include "common/type_aliases.hpp"
#include "options/options_block/misc_help_version_options_block.hpp"
#include "test/predicate/files_equal.hpp"
#include "test/predicate/string_matches_file.hpp"

#include <regex>
#include <sstream>

namespace cath { namespace test { } }

using namespace cath::clust;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::test;
using namespace std::literals::string_literals;

using boost::filesystem::path;
using boost::range::join;
using std::istringstream;
using std::ostringstream;
using std::regex;
using std::smatch;
using std::string;

namespace cath {
	namespace test {

		/// \brief The cluster_mapper_test_suite_fixture to assist in testing calc_hit_list
		struct cluster_mapper_test_suite_fixture : protected map_clusters_fixture {
		protected:
			~cluster_mapper_test_suite_fixture() noexcept {
				try {
					BOOST_CHECK_EQUAL( err_ss.str(), "" );
				}
				catch (...) {
				}
			}

			/// \brief Call perform_map_clusters() with the specified arguments preceded by a pseudo-program-name
			///        and the fixture's i/o streams.
			void execute_perform_map_clusters(const str_vec &prm_arguments ///< The arguments to pass to perform_map_clusters(), preceded by a pseudo-program-name
			                                  ) {
				// std::cerr << "Options are :\n\t\"" << boost::algorithm::join( prm_arguments, "\"\n\t\"" ) << "\n";

				const auto progname_rng = { "pseudo_program_name"s };
				perform_map_clusters(
					copy_build<str_vec>( join(
						progname_rng,
						prm_arguments
					) ),
					input_ss, output_ss, err_ss, parse_sources::CMND_LINE_ONLY
				);
			}

			/// \brief The input stream to use in the tests
			istringstream   input_ss;

			/// \brief The output stream to use in the tests
			ostringstream   output_ss;

			/// \brief The error stream to use in the tests
			ostringstream   err_ss;

			// /// \brief An output stream to which logging can be sent
			// ostringstream   log_ss;

			/// \brief A temporary temp_file
			const temp_file TEMP_TEST_FILE{ ".cath_cluster_mapper__test_file.%%%%-%%%%-%%%%-%%%%" };

			/// \brief The path of the TEMP_TEST_FILE
			const path      TEMP_TEST_FILE_FILENAME = get_filename( TEMP_TEST_FILE );
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(cath_cluster_mapper_test_suite, cluster_mapper_test_suite_fixture)

BOOST_AUTO_TEST_CASE(fails_if_given_no_options) {
	// Given an input stream containing the input
	input_ss.str( eg_input_str() );

	// When calling perform_map_clusters without options
	execute_perform_map_clusters( { } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(Must specify an input file)" } ) );
}

BOOST_AUTO_TEST_CASE(gives_correct_help_usage) {
	// When calling perform_map_clusters with options: help
	execute_perform_map_clusters( { "--" + misc_help_version_options_block::PO_HELP } );

	// Then expect the correct help message in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), help_usage_file() );
}

BOOST_AUTO_TEST_CASE(processes_from_stdin_to_stdout) {
	// Given an input stream containing the input
	input_ss.str( eg_input_str() );

	// When calling perform_map_clusters with options: - (a dash to read from the input stream)
	execute_perform_map_clusters( { "-" } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_renumber_only_result_file() );
}

BOOST_AUTO_TEST_CASE(processes_from_file_to_stdout) {
	// When calling perform_map_clusters with options: an input file
	execute_perform_map_clusters( { eg_input_file().string() } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_renumber_only_result_file() );
}

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

BOOST_AUTO_TEST_CASE(append_batch_id_works) {
	// When calling perform_map_clusters with options: an input file, the output file option and an output file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_output_options_block::PO_APPEND_BATCH_ID, "agifttome" } );

	// Then expect:
	//  * an empty output stream and
	//  * the output file containing the correct output
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_append_batch_id_result_file() );
}


BOOST_AUTO_TEST_SUITE(overlap_thresholds)

BOOST_AUTO_TEST_CASE(rejects_dom_overlap_option_if_not_mapping_from) {
	// When calling perform_map_clusters with options: an input file and a domain overlap
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL, "60" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(Cannot specify mapping threshold options)" } ) );
}

BOOST_AUTO_TEST_CASE(rejects_clust_overlap_option_if_not_mapping_from) {
	// When calling perform_map_clusters with options: an input file and a cluster overlap
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL, "60" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(Cannot specify mapping threshold options)" } ) );
}



BOOST_AUTO_TEST_SUITE(out_of_range)

BOOST_AUTO_TEST_CASE(rejects_dom_overlap_less_than_50) {
	// When calling perform_map_clusters with options: an input file and a too-small domain overlap
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL, "49.99" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(mapping fraction is out of range)" } ) );
}

BOOST_AUTO_TEST_CASE(rejects_clust_overlap_less_than_50) {
	// When calling perform_map_clusters with options: an input file and a too-small cluster overlap
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL, "49.99" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(mapping fraction is out of range)" } ) );
}

BOOST_AUTO_TEST_CASE(rejects_dom_overlap_more_than_100) {
	// When calling perform_map_clusters with options: an input file and a too-large domain overlap
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL, "100.1" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(mapping fraction is out of range)" } ) );
}

BOOST_AUTO_TEST_CASE(rejects_clust_overlap_more_than_100) {
	// When calling perform_map_clusters with options: an input file and a too-large cluster overlap
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL, "100.1" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(mapping fraction is out of range)" } ) );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(in_range)

BOOST_AUTO_TEST_CASE(accepts_dom_overlap_of_50) {
	// When calling perform_map_clusters with options: an input file, a map-from file and a domain overlap of 50%
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_mapfrom_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL, "50" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_mapfrom_dom_ol_50_result_file() );
}

BOOST_AUTO_TEST_CASE(accepts_clust_overlap_of_50) {
	// When calling perform_map_clusters with options: an input file, a map-from file and a cluster overlap of 50%
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_mapfrom_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL, "50" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_mapfrom_clust_ol_50_result_file() );
}

BOOST_AUTO_TEST_CASE(accepts_dom_overlap_of_100) {
	// When calling perform_map_clusters with options: an input file, a map-from file and a domain overlap of 100%
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_mapfrom_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL, "100" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_mapfrom_dom_ol_100_result_file() );
}

BOOST_AUTO_TEST_CASE(accepts_clust_overlap_of_100) {
	// When calling perform_map_clusters with options: an input file, a map-from file and a cluster overlap of 100%
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_mapfrom_file().string(),
		"--" + clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL, "100" } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_mapfrom_clust_ol_100_result_file() );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE_END()




BOOST_AUTO_TEST_CASE(fails_if_batch_id_when_using_batches) {
	// When calling perform_map_clusters with options: an input file, the append-batch-id flag and the batches flag
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_output_options_block::PO_APPEND_BATCH_ID, "agifttome",
		"--" + clustmap_input_options_block::PO_READ_BATCHES_FROM_INPUT } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(Cannot specify a batch ID for appending.*when reading batches from input)" } ) );
}


BOOST_AUTO_TEST_CASE(fails_if_map_from_file_when_using_batches) {
	// When calling perform_map_clusters with options: an input file, a map-from file and the batches flag
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_mapfrom_file().string(),
		"--" + clustmap_input_options_block::PO_READ_BATCHES_FROM_INPUT } );

	// Then expect the correct error message in the output stream
	BOOST_CHECK( regex_search( output_ss.str(), regex{ R"(Cannot specify a map-from cluster-membership file.*when reading batches from input)" } ) );
}


BOOST_AUTO_TEST_CASE(accepts_non_numeric_cluster_names_in_map_from) {
	// When calling perform_map_clusters with options: an input file, a map-from file that includes non-numeric cluster names
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_non_numeric_file().string()} );

	// Then expect the correct results output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_input_non_numeric_fromresult() );
}


BOOST_AUTO_TEST_CASE(accepts_non_numeric_cluster_names_in_map_to) {
	// When calling perform_map_clusters with options: an input file that includes non-numeric cluster names
	execute_perform_map_clusters( { eg_input_non_numeric_file().string() } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_input_non_numeric_toresult() );
}


/// \TODO Prefer a case that has some mapping results below 100% to better check the percentile stats
BOOST_AUTO_TEST_CASE(generates_correct_summary_file_when_mapping) {
	// When calling perform_map_clusters with options: an input file, a map-from file and a file to write a summary to
	execute_perform_map_clusters( { eg_input_non_numeric_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_mapfrom_file().string(),
		"--" + clustmap_output_options_block::PO_SUMMARISE_TO_FILE,      TEMP_TEST_FILE_FILENAME.string(),
		// "--" + clust_mapping_options_block::PO_MIN_EQUIV_DOM_OL,         "99",
		// "--" + clust_mapping_options_block::PO_MIN_EQUIV_CLUST_OL,       "99",
	} );

	// Then expect the correct output in the summary file (TEMP_TEST_FILE_FILENAME)
	BOOST_CHECK_FILES_EQUAL( TEMP_TEST_FILE_FILENAME, eg_summary_mapping_file() );
}

BOOST_AUTO_TEST_CASE(generates_correct_summary_file_when_renumbering) {
	// When calling perform_map_clusters with options: an input file and a file to write a summary to
	execute_perform_map_clusters( { eg_input_non_numeric_file().string(),
		"--" + clustmap_output_options_block::PO_SUMMARISE_TO_FILE, TEMP_TEST_FILE_FILENAME.string() } );

	// Then expect the correct output in the summary file (TEMP_TEST_FILE_FILENAME)
	BOOST_CHECK_FILES_EQUAL( TEMP_TEST_FILE_FILENAME, eg_renumbering_summary_file() );
}

// Check that the mapping correctly require that:
//  * > 0.2 of the candidate new cluster's entries map to the candidate old cluster and
//  * > 0.5 of the candidate new cluster's entries that map to some old cluster map to the candidate old cluster
//
// The correct answer is:
//  * to_cluster_w *should not* map to from_cluster_a (and to nothing else)
//  * to_cluster_x *should not* map to from_cluster_b (and to nothing else)
//  * to_cluster_y *should not* map to from_cluster_c (and to nothing else)
//  * to_cluster_z *should*     map to from_cluster_d (and to nothing else)
//
// The breakdown:
//
//   from_cluster_name | from_size | from_mapped | to_cluster_name | to_size | to_mapped || pair_mapped | frac_to | frac_to_mapped | decision
//  -------------------+-----------+-------------+-----------------+---------+-----------++-------------+---------+----------------+----------
//   from_cluster_a    |         2 |           2 | to_cluster_w    |      10 |         4 ||           2 |   0.2 ✗ |         0.5  ✗ |    ✗
//   from_cluster_b    |         5 |           5 | to_cluster_x    |      10 |        10 ||           5 |   0.5 ✓ |         0.5  ✗ |    ✗
//   from_cluster_c    |         2 |           2 | to_cluster_y    |      10 |         2 ||           2 |   0.2 ✗ |         1.0  ✓ |    ✗
//   from_cluster_d    |         3 |           3 | to_cluster_z    |      10 |         4 ||           3 |   0.3 ✓ |         0.75 ✓ |    ✓
//
BOOST_AUTO_TEST_CASE(handles_to_clust_ol_thresholds_correctly) {
	// When calling perform_map_clusters with options: an input file and a map-from file with data to test to-cluster thresholds
	execute_perform_map_clusters( { to_clust_ol_thresholds_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, to_clust_ol_thresholds_mapfrom_file().string(),
	} );

	// Then expect the correct results in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), to_clust_ol_thresholds_result_file() );
}


BOOST_AUTO_TEST_SUITE(clashes_duplicates_etc)

BOOST_AUTO_TEST_CASE(handles_clashing_segments_in_input) {
	// When calling perform_map_clusters with options: a clashing_segments input file
	execute_perform_map_clusters( { eg_input_clashing_segments_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(clashes with a previous entry)" } ), 2 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_clashing_segments_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a clashing_segments map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_clashing_segments_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(clashes with a previous entry)" } ), 2 );
	err_ss.str( "" );
}




BOOST_AUTO_TEST_CASE(handles_clashing_segments_w_diff_names_in_input) {
	// When calling perform_map_clusters with options: a clashing_segments_w_diff_names input file
	execute_perform_map_clusters( { eg_input_clashing_segments_w_diff_names_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK( regex_search( err_ss.str(), regex{ R"(clashes with a previous entry)" } ) );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_clashing_segments_w_diff_names_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a clashing_segments_w_diff_names map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_clashing_segments_w_diff_names_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK( regex_search( err_ss.str(), regex{ R"(clashes with a previous entry)" } ) );
	err_ss.str( "" );
}



BOOST_AUTO_TEST_CASE(handles_mixed_wcds_and_segments_in_input) {
	// When calling perform_map_clusters with options: a mixed_wcds_and_segments input file
	execute_perform_map_clusters( { eg_input_mixed_wcds_and_segments_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(clashes with a previous entry)" } ), 1 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_mixed_wcds_and_segments_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a mixed_wcds_and_segments map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_mixed_wcds_and_segments_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( err_ss.str(), "" ); // Ideally, it'd be good if this detected the clash across clusters but
	                                       // at present the map-from data is stored separately by cluster so
	                                       // it doesn't seem worth the additional effort
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_eg_input_mixed_wcds_and_segments_within_cluster_file_in_input) {
	// When calling perform_map_clusters with options: a mixed_wcds_and_segments_within_cluster input file
	execute_perform_map_clusters( { eg_input_mixed_wcds_and_segments_within_cluster_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(clashes with a previous entry)" } ), 1 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_eg_input_mixed_wcds_and_segments_within_cluster_file_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a mixed_wcds_and_segments_within_cluster map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_mixed_wcds_and_segments_within_cluster_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(clashes with a previous entry)" } ), 1 );
	err_ss.str( "" );
}




BOOST_AUTO_TEST_CASE(handles_repeated_segments_in_input) {
	// When calling perform_map_clusters with options: a repeated_segments input file
	execute_perform_map_clusters( { eg_input_repeated_segments_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(duplicates a previous entry)" } ), 1 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_repeated_segments_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a repeated_segments map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_repeated_segments_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(duplicates a previous entry)" } ), 1 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_repeated_segments_w_diff_names_in_input) {
	// When calling perform_map_clusters with options: a repeated_segments_w_diff_names input file
	execute_perform_map_clusters( { eg_input_repeated_segments_w_diff_names_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(duplicates a previous entry)" } ), 1 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_repeated_segments_w_diff_names_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a repeated_segments_w_diff_names map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_repeated_segments_w_diff_names_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(duplicates a previous entry)" } ), 1 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_repeated_wcds_in_input) {
	// When calling perform_map_clusters with options: a repeated_wcds input file
	execute_perform_map_clusters( { eg_input_repeated_wcds_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(duplicates a previous entry)" } ), 1 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_repeated_wcds_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a repeated_wcds map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_repeated_wcds_file().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_NE( output_ss.str(), "" );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(duplicates a previous entry)" } ), 1 );
	err_ss.str( "" );
}





BOOST_AUTO_TEST_CASE(handles_start_at_zero_in_input) {
	// When calling perform_map_clusters with options: a zero_start input file
	execute_perform_map_clusters( { eg_input_zero_start().string() } );

	// Then expect the correct output
	BOOST_CHECK_EQUAL( output_ss.str(), "# cluster-id suggested-name\nclust_a 1\n" );
}

BOOST_AUTO_TEST_CASE(handles_start_at_zero_in_input_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a zero_start map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_zero_start().string() } );

	// Then expect the correct output
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_simple_mapfrom_result_file() );
}




BOOST_AUTO_TEST_CASE(handles_backward_segment_in_input) {
	// When calling perform_map_clusters with options: a backward segment input file
	execute_perform_map_clusters( { eg_input_backward_segment().string()  } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), empty_result_file() );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(parsing segments from entry.*start residue before the stop residue)" } ), 1 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_backward_segment_in_input_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a backward segment map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_backward_segment().string()  } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_simple_mapfrom_result_file() );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(parsing segments from entry.*start residue before the stop residue)" } ), 1 );
	err_ss.str( "" );
}



BOOST_AUTO_TEST_CASE(handles_misordered_segments_in_input) {
	// When calling perform_map_clusters with options: a misordered segments input file
	execute_perform_map_clusters( { eg_input_misordered_segments().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), empty_result_file() );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(parsing segments from entry.*building segments from bounds, found preceding stop after start)" } ), 1 );
	err_ss.str( "" );
}

BOOST_AUTO_TEST_CASE(handles_misordered_segments_in_input_in_map_from) {
	// When calling perform_map_clusters with options: an input file and a misordered segments map-from file
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_misordered_segments().string() } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_simple_mapfrom_result_file() );
	BOOST_CHECK_EQUAL( regex_count( err_ss.str(), regex{ R"(parsing segments from entry.*building segments from bounds, found preceding stop after start)" } ), 1 );
	err_ss.str( "" );
}






BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_CASE(handles_batch) {
	// When calling perform_map_clusters with options: an batch input file and the --read-batches-from-input flag
	execute_perform_map_clusters( { eg_batch_input_file().string(),
		"--" + clustmap_input_options_block::PO_READ_BATCHES_FROM_INPUT } );

	// Then expect the correct error message in the error stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_batch_result_file() );
}


BOOST_AUTO_TEST_CASE(provides_entry_level_output) {
	// When calling perform_map_clusters with options: an input file, a map-from file and the --print-entry-results flag
	execute_perform_map_clusters( { eg_input_file().string(),
		"--" + clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE, eg_input_mapfrom_file().string(),
		"--" + clustmap_output_options_block::PO_PRINT_DOMAIN_MAPPING } );

	// Then expect the correct results in the output stream
	BOOST_CHECK_STRING_MATCHES_FILE( output_ss.str(), eg_entry_mapping_result_file() );
}


BOOST_AUTO_TEST_SUITE_END()
