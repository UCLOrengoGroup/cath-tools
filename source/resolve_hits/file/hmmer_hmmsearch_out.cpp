/// \file
/// \brief The hmmer_hmmsearch_out definitions

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

#include "hmmer_hmmsearch_out.hpp"

#include "common/file/open_fstream.hpp"
#include "resolve_hits/file/detail/hmmsearch_parser.hpp"
#include "resolve_hits/read_and_process_hits/read_and_process_mgr.hpp"

#include <fstream>

using namespace cath;
using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;
using namespace cath::seq;

using boost::filesystem::path;
using std::ifstream;
using std::istream;

/// \brief Parse hmmsearch output data from the specified file and pass the hits to the specified read_and_process_mgr
void cath::rslv::parse_hmmsearch_out_file(read_and_process_mgr &arg_read_and_process_mgr, ///< The read_and_process_mgr to which the hits should be passed for processing
                                          const path           &arg_hmmsearch_out_file,   ///< The file from which the HMMER domain hits table data should be parsed
                                          const bool           &arg_apply_cath_policies,  ///< Whether to apply CATH-specific policies
                                          const residx_t       &arg_min_gap_length,       ///< The minimum length that an alignment gap can have to be considered a gap
                                          const bool           &arg_output_hmmsearch_aln  ///< Whether to parse/output hmmsearch output alignment information
                                          ) {
	ifstream the_ifstream;
	open_ifstream( the_ifstream, arg_hmmsearch_out_file );

	parse_hmmsearch_out(
		arg_read_and_process_mgr,
		the_ifstream,
		arg_apply_cath_policies,
		arg_min_gap_length,
		arg_output_hmmsearch_aln
	);

	the_ifstream.close();
}

/// \brief Parse hmmsearch output data from the specified input stream and pass the hits to the specified read_and_process_mgr
void cath::rslv::parse_hmmsearch_out(read_and_process_mgr &arg_read_and_process_mgr, ///< The read_and_process_mgr to which the hits should be passed for processing
                                     istream              &arg_input_stream,         ///< The istream from which the HMMER domain hits table data should be parsed
                                     const bool           &arg_apply_cath_policies,  ///< Whether to apply CATH-specific policies
                                     const residx_t       &arg_min_gap_length,       ///< The minimum length that an alignment gap can have to be considered a gap
                                     const bool           &arg_parse_hmmsearch_aln   ///< Whether to parse/output hmmsearch output alignment information
                                     ) {
	hmmsearch_parser parser{ arg_input_stream };

	arg_read_and_process_mgr.process_all_outstanding();

	// Store the query IDs seen so far if the crh_filter_spec specifies a limit on the number of queries
	query_id_recorder seen_query_ids;

	while ( ! parser.end_of_istream() ) {
		parser.advance_line_to_block();

		if ( parser.end_of_istream() ) {
			break;
		}

		parser.advance_line();

		// If this query ID should be skipped, then skip this entry.
		// The function also updates seen_query_ids if not skipping this query ID
		if ( should_skip_query_and_update( arg_read_and_process_mgr, parser.get_query_id(), seen_query_ids ) ) {
			continue;
		}

		if ( parser.line_is_summary_header() ) {
			parser.parse_summary_from_header(
				arg_apply_cath_policies,
				arg_read_and_process_mgr.get_filter_spec()
			);
			parser.advance_line_until_next_aln();

			while ( ! parser.line_is_at_block() && ! parser.line_is_at_pipeline_stats() ) {
				if ( parser.line_is_at_aln() ) {
					if ( ! parser.alignment_is_empty() ) {
						parser.finish_alignment(
							arg_read_and_process_mgr,
							arg_apply_cath_policies,
							arg_min_gap_length,
							arg_parse_hmmsearch_aln
						);
					}
					parser.advance_line();
				}
				parser.parse_alignment_section();
			}
			if ( ! parser.alignment_is_empty() ) {
				parser.finish_alignment(
					arg_read_and_process_mgr,
					arg_apply_cath_policies,
					arg_min_gap_length,
					arg_parse_hmmsearch_aln
				);
			}
		}
	}

	arg_read_and_process_mgr.process_all_outstanding();
}
