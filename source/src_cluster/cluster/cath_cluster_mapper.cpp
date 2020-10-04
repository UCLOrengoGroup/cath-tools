/// \file
/// \brief The cath_cluster_mapper class definitions

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

#include "cath_cluster_mapper.hpp"

#include <boost/filesystem/path.hpp>

#include "cluster/clustmap_options.hpp"
#include "cluster/detail/mapping_job.hpp"
#include "cluster/file/cluster_membership_file.hpp"
#include "cluster/map/aggregate_map_results.hpp"
#include "cluster/map/map_clusters.hpp"
#include "cluster/map/map_results.hpp"
#include "cluster/new_cluster_data.hpp"
#include "cluster/old_cluster_data.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/container/id_of_str_bidirnl.hpp"
#include "common/file/ofstream_list.hpp"
#include "common/file/open_fstream.hpp"
#include "common/file/path_or_istream.hpp"
#include "common/logger.hpp"
#include "common/optional/make_optional_if.hpp"
#include "options/executable/executable_options.hpp"

#include <fstream>
#include <functional>

using namespace cath::clust::detail;
using namespace cath::common;
using namespace cath::opts;

using boost::filesystem::path;
using std::ifstream;
using std::istream;
using std::ostream;

/// \brief Perform map-clusters according to the specified arguments strings with the specified i/o streams
void cath::clust::perform_map_clusters(const str_vec       &args,             ///< The arguments strings specifying the map-clusters action to perform
                                       istream             &prm_istream,      ///< The input stream
                                       ostream             &prm_stdout,       ///< The output stream
                                       ostream             &prm_stderr,       ///< The error stream
                                       const parse_sources &prm_parse_sources ///< The sources from which options should be parsed
                                       ) {
	perform_map_clusters(
		make_and_parse_options<clustmap_options>( args, prm_parse_sources ),
		prm_istream,
		prm_stdout,
		prm_stderr
	);
}

/// \brief Perform map-clusters according to the specified clustmap_options with the specified i/o streams
void cath::clust::perform_map_clusters(const clustmap_options &prm_opts,    ///< The clustmap_options specifying the map-clusters action to perform
                                       istream                &prm_istream, ///< The input stream
                                       ostream                &prm_stdout,  ///< The output stream
                                       ostream                &prm_stderr   ///< The error stream
                                       ) {
	// If the options are invalid or specify to do_nothing, then just return
	const auto &error_or_help_string = prm_opts.get_error_or_help_string();
	if ( error_or_help_string ) {
		prm_stdout << *error_or_help_string;
		return;
	}

	perform_map_clusters(
		prm_opts.get_clustmap_input_spec(),
		prm_opts.get_clust_mapping_spec(),
		prm_opts.get_clustmap_output_spec(),
		prm_istream,
		prm_stdout,
		prm_stderr
	);
}

/// \brief Perform map-clusters according to the specified crh_spec with the specified i/o streams
void cath::clust::perform_map_clusters(const clustmap_input_spec   &prm_input_spec,   ///< The clustmap_input_spec specifying the map-clusters action to perform
                                       const clust_mapping_spec    &prm_mapping_spec, ///< The clustmap_mapping_spec specifying the map-clusters action to perform
                                       const clustmap_output_spec  &prm_output_spec,  ///< The clustmap_output_spec specifying the map-clusters action to perform
                                       istream                     &prm_istream,      ///< The input stream
                                       ostream                     &prm_stdout,       ///< The output stream
                                       ostream                     &prm_stderr        ///< The error stream
                                       ) {
	const path     &working_clustmemb_file  = prm_input_spec.get_working_clustmemb_file();
	const bool     &read_batches_from_input = prm_input_spec.get_read_batches_from_input();

	// Organise the input stream
	if ( working_clustmemb_file != "-" && ! exists( working_clustmemb_file ) ) {
		logger::log_and_exit(
			logger::return_code::NO_SUCH_FILE,
			"No such map-clusters input data file \"" + working_clustmemb_file.string() + "\""
		);
	}

	path_or_istream istream_wrapper{ prm_istream };
	istream_wrapper.set_path( working_clustmemb_file ).get_istream();

	const mapping_job_vec jobs = read_batches_from_input ? read_batch_mapping_file(
	                                                     	istream_wrapper.set_path( working_clustmemb_file ).get_istream()
	                                                     )
	                                                     : mapping_job_vec{ {
	                                                     	mapping_job{
	                                                     		prm_output_spec.get_append_batch_id(),
	                                                     		working_clustmemb_file,
	                                                     		prm_input_spec.get_map_from_clustmemb_file()
	                                                     	}
	                                                     } };
	istream_wrapper.close();

	aggregate_map_results agg_results{ prm_mapping_spec };

	ofstream_list out_list{ prm_stdout };
	auto the_ostreams = out_list.open_ofstreams(
		path_vec{ { prm_output_spec.get_output_to_file().value_or( out_list.get_flag() ) } }
	);

	// \TODO Come C++17 and structured bindings, use here
	for (const size_t &job_idx : indices( jobs.size() ) ) {
		const mapping_job &job                    = jobs[ job_idx ];
		const auto        &job_batch_id           = job.get_batch_id();
		const auto        &job_new_clustmemb_file = job.get_new_cluster_membership_file();
		const auto        &job_old_clustmemb_file = job.get_old_cluster_membership_file();

		id_of_str_bidirnl seq_ider;

		auto &the_istream = istream_wrapper.set_path( job_new_clustmemb_file ).get_istream();
		const new_cluster_data new_to_clusters = parse_new_membership( the_istream, seq_ider, ref( prm_stderr ) );
		istream_wrapper.close();

		const old_cluster_data_opt old_from_clusters = make_optional_if_fn(
			static_cast<bool>( job_old_clustmemb_file ),
			[&] { return parse_old_membership( *job_old_clustmemb_file, seq_ider, ref( prm_stderr )  ); }
		);

		const auto results = map_clusters(
			old_from_clusters,
			new_to_clusters,
			prm_mapping_spec,
			make_optional_if( prm_output_spec.get_print_domain_mapping(), the_ostreams.front() )
		);

		if ( ! prm_output_spec.get_print_domain_mapping() ) {
			for (auto &the_ostream : the_ostreams) {
				the_ostream.get() << results_string( old_from_clusters, new_to_clusters, results, job_batch_id, ( job_idx == 0 ) );
				// the_ostream.get() << longer_results_string( old_from_clusters, new_to_clusters, results, job_batch_id );
			}
		}

		if ( old_from_clusters ) {
			agg_results.add_map_results( results, *old_from_clusters, new_to_clusters );
		}
	}

	if ( prm_output_spec.get_summarise_to_file() ) {
		write_markdown_summary_string_to_file( *prm_output_spec.get_summarise_to_file(), agg_results );
	}
}
