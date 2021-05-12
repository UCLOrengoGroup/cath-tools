/// \file
/// \brief The parse_hmmer_out definitions

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

#include "parse_hmmer_out.hpp"

#include <filesystem>
#include <fstream>
#include <memory>

#include "cath/common/file/open_fstream.hpp"
#include "cath/resolve_hits/file/detail/hmmer_parser.hpp"
#include "cath/resolve_hits/read_and_process_hits/hits_processor/gather_hits_processor.hpp"
#include "cath/resolve_hits/read_and_process_hits/hits_processor/hits_processor_list.hpp"
#include "cath/resolve_hits/read_and_process_hits/read_and_process_mgr.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::rslv;
using namespace ::cath::rslv::detail;
using namespace ::cath::seq;

using ::std::filesystem::path;
using ::std::ifstream;
using ::std::istream;

/// \brief Parse HMMER output data from the specified file and pass the hits to the specified read_and_process_mgr
void cath::rslv::parse_hmmer_out_file(read_and_process_mgr &prm_read_and_process_mgr, ///< The read_and_process_mgr to which the hits should be passed for processing
                                      const path           &prm_hmmer_out_file,       ///< The file from which the HMMER domain hits table data should be parsed
                                      const hmmer_format   &prm_hmmer_format,         ///< The HMMER format to parse
                                      const bool           &prm_apply_cath_policies,  ///< Whether to apply CATH-specific policies
                                      const residx_t       &prm_min_gap_length,       ///< The minimum length that an alignment gap can have to be considered a gap
                                      const bool           &prm_output_hmmer_aln      ///< Whether to parse/output HMMER output alignment information
                                      ) {
	ifstream the_ifstream = open_ifstream( prm_hmmer_out_file );

	parse_hmmer_out(
		prm_read_and_process_mgr,
		the_ifstream,
		prm_hmmer_format,
		prm_apply_cath_policies,
		prm_min_gap_length,
		prm_output_hmmer_aln
	);

	the_ifstream.close();
}

/// \brief Parse HMMER output data from the specified input stream and pass the hits to the specified read_and_process_mgr
void cath::rslv::parse_hmmer_out(read_and_process_mgr &prm_read_and_process_mgr, ///< The read_and_process_mgr to which the hits should be passed for processing
                                 istream              &prm_input_stream,         ///< The istream from which the HMMER domain hits table data should be parsed
                                 const hmmer_format   &prm_hmmer_format,         ///< The HMMER format to parse
                                 const bool           &prm_apply_cath_policies,  ///< Whether to apply CATH-specific policies
                                 const residx_t       &prm_min_gap_length,       ///< The minimum length that an alignment gap can have to be considered a gap
                                 const bool           &prm_parse_hmmer_aln       ///< Whether to parse/output HMMER output alignment information
                                 ) {
	hmmer_parser parser{ prm_hmmer_format, prm_input_stream };

	prm_read_and_process_mgr.process_all_outstanding();

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
		if ( should_skip_query_and_update( prm_read_and_process_mgr, parser.get_query_id(), seen_query_ids ) ) {
			continue;
		}

		if ( parser.line_is_summary_header() ) {
			parser.parse_summary_from_header(
				prm_apply_cath_policies,
				prm_read_and_process_mgr.get_filter_spec()
			);
			parser.advance_line_until_next_aln();

			while ( ! parser.line_is_at_block() && ! parser.line_is_at_pipeline_stats() ) {
				if ( parser.line_is_at_aln() ) {
					if ( ! parser.alignment_is_empty() ) {
						parser.finish_alignment(
							prm_read_and_process_mgr,
							prm_apply_cath_policies,
							prm_min_gap_length,
							prm_parse_hmmer_aln
						);
					}
					parser.advance_line();
				}
				parser.parse_alignment_section();
			}
			if ( ! parser.alignment_is_empty() ) {
				parser.finish_alignment(
					prm_read_and_process_mgr,
					prm_apply_cath_policies,
					prm_min_gap_length,
					prm_parse_hmmer_aln
				);
			}
		}
	}

	prm_read_and_process_mgr.process_all_outstanding();
}

/// \brief Parse HMMER output data from the specified file and pass the hits to the specified read_and_process_mgr
str_calc_hit_list_pair_vec cath::rslv::parse_hmmer_out_file(const path         &prm_hmmer_out_file,      ///< The file from which the HMMER domain hits table data should be parsed
                                                            const hmmer_format &prm_hmmer_format,        ///< The HMMER format to parse
                                                            const bool         &prm_apply_cath_policies, ///< Whether to apply CATH-specific policies
                                                            const residx_t     &prm_min_gap_length,      ///< The minimum length that an alignment gap can have to be considered a gap
                                                            const bool         &prm_output_hmmer_aln     ///< Whether to parse/output HMMER output alignment information
                                                            ) {
	// gather_hits_processor the_processor;
	str_calc_hit_list_pair_vec results;
	hits_processor_list hits_processors;
	hits_processors.add_processor( gather_hits_processor{ results } );
	read_and_process_mgr the_read_and_process_mgr{
		hits_processors,
		crh_filter_spec{}
	};
	parse_hmmer_out_file(
		the_read_and_process_mgr,
		prm_hmmer_out_file,
		prm_hmmer_format,
		prm_apply_cath_policies,
		prm_min_gap_length,
		prm_output_hmmer_aln
	);
	return results;
}
