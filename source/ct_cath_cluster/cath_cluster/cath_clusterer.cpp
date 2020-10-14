/// \file
/// \brief The cath_clusterer class definitions

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

#include "cath_clusterer.hpp"

#include <boost/log/trivial.hpp>

#include "cath_cluster/options/cath_cluster_options.hpp"
#include "clustagglom/calc_complete_linkage_merge_list.hpp"
#include "clustagglom/file/dissimilarities_file.hpp"
#include "clustagglom/file/names_file.hpp"
#include "clustagglom/get_sorting_scores.hpp"
#include "clustagglom/hierarchy.hpp"
#include "clustagglom/link_dirn.hpp"
#include "clustagglom/links.hpp"
#include "clustagglom/make_clusters_from_merges.hpp"
#include "clustagglom/merge.hpp"
#include "common/container/id_of_str_bidirnl.hpp"
#include "common/file/ofstream_list.hpp"
#include "common/file/open_fstream.hpp"
#include "common/file/path_or_istream.hpp"
#include "common/logger.hpp"

#include <fstream>

using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::boost::filesystem::path;
using ::boost::log::trivial::warning;
using ::boost::make_optional;
using ::std::ifstream;
using ::std::istream;
using ::std::move;
using ::std::ostream;

/// \brief Perform clustering according to the specified arguments strings with the specified i/o streams
void cath::clust::perform_cluster(const str_vec       &args,             ///< The arguments strings specifying the clustering action to perform
                                  istream             &prm_istream,      ///< The input stream
                                  ostream             &prm_stdout,       ///< The output stream
                                  const parse_sources &prm_parse_sources ///< The sources from which options should be parsed
                                  ) {
	perform_cluster(
		make_and_parse_options<cath_cluster_options>( args, prm_parse_sources ),
		prm_istream,
		prm_stdout
	);
}

/// \brief Perform clustering according to the specified cath_cluster_options with the specified i/o streams
void cath::clust::perform_cluster(const cath_cluster_options &prm_opts,    ///< The cath_cluster_options specifying the clustering action to perform
                                  istream                    &prm_istream, ///< The input stream
                                  ostream                    &prm_stdout   ///< The output stream
                                  ) {
	// If the options are invalid or specify to do_nothing, then just return
	const auto &error_or_help_string = prm_opts.get_error_or_help_string();
	if ( error_or_help_string ) {
		prm_stdout << *error_or_help_string;
		return;
	}

	const auto         &in_spec        = prm_opts.get_cath_cluster_input_spec();
	const auto         &clust_spec     = prm_opts.get_cath_cluster_clustering_spec();
	const auto         &out_spec       = prm_opts.get_cath_cluster_output_spec();
	const link_dirn    &the_link_dirn  = in_spec.get_link_dirn();
	const size_t       &column_idx     = in_spec.get_column_idx();
	const path_opt     &links_infile   = in_spec.get_links_infile();
	const str_opt       level_warning  = get_dissim_sort_warning( clust_spec, the_link_dirn );
	const strength_vec  cutoffs        = get_sorted_dissims     ( clust_spec, the_link_dirn );
	const strength      the_max_dissim = get_max_dissim         ( clust_spec, the_link_dirn );

	if ( level_warning ) {
		BOOST_LOG_TRIVIAL( warning ) << *level_warning;
	}

	// Organise the input stream
	path_or_istream istream_wrapper{ prm_istream };
	if ( file_is_missing( istream_wrapper, *links_infile ) ) {
		logger::log_and_exit(
			logger::return_code::NO_SUCH_FILE,
			"No such links input data file \"" + links_infile->string() + "\""
		);
	}
	istream &the_istream = istream_wrapper.set_path( *links_infile ).get_istream();

	id_of_str_bidirnl the_name_ider;
	const bool has_names_file = static_cast<bool>( in_spec.get_names_infile() );

	// If there is a names file, want to parse that before parsing the links
	// but if there isn't a links file, want to parse the links before getting the sorting scores
	const doub_vec props           = has_names_file ? parse_names( *in_spec.get_names_infile(), the_name_ider ) : doub_vec{};
	const auto     dissims         = parse_dissimilarities( the_istream, the_name_ider, the_link_dirn, column_idx );
	const size_vec sorting_indices = [&] {
		if ( has_names_file ) {
			return get_sorting_scores( the_name_ider, props );
		}
		BOOST_LOG_TRIVIAL( warning )
			<< "No names file has been specified. You are recommended to specify a names file ("
			<< cath_cluster_input_options_block::PO_NAMES_INFILE
			<< ") to ensure singletons don't get missed.";
		return get_sorting_scores( the_name_ider );
	} ();

	// In principle, the dissims can be made non-const and then passed to
	// calc_complete_linkage_merge_list() via a std::move() iff it isn't going to
	// be required later on for write_spanning_trees(). In practice, you can't
	// conditionally move() based on a run-time condition.
	//
	// But another way of avoiding the copy when not writing a spanning tree is
	// to always move here but conditionally take a copy for use in write_spanning_trees()
	// beforehand.
	const auto     merges          = calc_complete_linkage_merge_list(
		dissims,
		sorting_indices,
		the_max_dissim
	);

	ofstream_list ofstreams{ prm_stdout };

	// If merges output has been requested then write it
	if ( out_spec.get_merges_to_file() ) {
		write_merge_list( open_ofstream( ofstreams, *out_spec.get_merges_to_file() ).get(), merges );
	}

	// Build a hierarchy from the merges
	// (Should be O(n) so substantially cheaper than forming the merges)
	const hierarchy the_hierarchy = make_clusters_from_merges_and_sort(
		merges,
		sorting_indices,
		cutoffs
	);

	// If clusters output has been requested then write it
	const path_opt clusters_to_file = ( get_num_output_paths( out_spec ) == 0 )
	                                  ? make_optional( ofstreams.get_flag() )
	                                  : out_spec.get_clusters_to_file();
	if ( clusters_to_file ) {
		write_cluster(
			open_ofstream( ofstreams, *clusters_to_file ).get(),
			the_hierarchy,
			the_name_ider
		);
	}

	// If spanning trees has been requested then write it
	if ( out_spec.get_clust_spans_to_file() ) {
		write_spanning_trees(
			open_ofstream( ofstreams, *out_spec.get_clust_spans_to_file() ).get(),
			the_hierarchy,
			the_name_ider,
			dissims
		);
	}
	// If reps output has been requested then write it
	if ( out_spec.get_reps_to_file() ) {
		write_reps(
			open_ofstream( ofstreams, *out_spec.get_reps_to_file() ).get(),
			the_hierarchy,
			the_name_ider
		);
	}
	// If sorted links output has been requested then write it
	if ( out_spec.get_sorted_links_to_file() ) {
		write_ordered_links(
			open_ofstream( ofstreams, *out_spec.get_sorted_links_to_file() ).get(),
			dissims,
			the_name_ider,
			sorting_indices
		);
	}
}
