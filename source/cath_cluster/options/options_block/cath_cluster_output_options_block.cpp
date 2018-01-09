/// \file
/// \brief The cath_cluster_output_options_block class definitions

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

#include "cath_cluster_output_options_block.hpp"

#include <boost/filesystem/path.hpp>

#include "common/clone/make_uptr_clone.hpp"

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::boost::filesystem::path;
using ::boost::none;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::string;
using ::std::unique_ptr;

/// \brief The option name for an optional file to which clusters should be written
const string cath_cluster_output_options_block::PO_CLUSTERS_TO_FILE     { "clusters-to-file"     };

/// \brief The option name for an optional file to which merges should be written
const string cath_cluster_output_options_block::PO_MERGES_TO_FILE       { "merges-to-file"       };

/// \brief The option name for an optional file to which clust_spans should be written
const string cath_cluster_output_options_block::PO_CLUST_SPANS_TO_FILE  { "clust-spans-to-file"  };

/// \brief The option name for an optional file to which reps should be written
const string cath_cluster_output_options_block::PO_REPS_TO_FILE         { "reps-to-file"         };

/// \brief The option name for an optional file to which sorted_links should be written
const string cath_cluster_output_options_block::PO_SORTED_LINKS_TO_FILE { "sorted-links-to-file" };

/// \brief A standard do_clone method
unique_ptr<options_block> cath_cluster_output_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string cath_cluster_output_options_block::do_get_block_name() const {
	return "Output";
}

/// \brief Add this block's options to the provided options_description
void cath_cluster_output_options_block::do_add_visible_options_to_description(options_description &arg_desc,           ///< The options_description to which the options are added
                                                                              const size_t        &/*arg_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                              ) {
	const string file_varname   { "<file>" };

	const auto clusters_to_file_notifier     = [&] (const path &x) { the_spec.set_clusters_to_file    ( x ); };
	const auto merges_to_file_notifier       = [&] (const path &x) { the_spec.set_merges_to_file      ( x ); };
	const auto clust_spans_to_file_notifier  = [&] (const path &x) { the_spec.set_clust_spans_to_file ( x ); };
	const auto reps_to_file_notifier         = [&] (const path &x) { the_spec.set_reps_to_file        ( x ); };

	arg_desc.add_options()
		(
			PO_CLUSTERS_TO_FILE.c_str(),
			value<path>()
				->value_name   ( file_varname                 )
				->notifier     ( clusters_to_file_notifier    ),
			( "Write the clustering to file "
			  + file_varname
			  + " (or '-' for stdout)" ).c_str()
		)
		(
			PO_MERGES_TO_FILE.c_str(),
			value<path>()
				->value_name   ( file_varname                 )
				->notifier     ( merges_to_file_notifier      ),
			( "Write the ordered list of merges to file "
			  + file_varname
			  + " (or '-' for stdout)" ).c_str()
		)
		(
			PO_CLUST_SPANS_TO_FILE.c_str(),
			value<path>()
				->value_name   ( file_varname                 )
				->notifier     ( clust_spans_to_file_notifier ),
			( "Write links that form spanning trees for each cluster to file "
			  + file_varname
			  + " (or '-' for stdout)" ).c_str()
		)
		(
			PO_REPS_TO_FILE.c_str(),
			value<path>()
				->value_name   ( file_varname                 )
				->notifier     ( reps_to_file_notifier        ),
			( "Write the list of representatives to file "
			  + file_varname
			  + " (or '-' for stdout)" ).c_str()
		);
}

/// \brief Add any hidden options to the provided options_description
void cath_cluster_output_options_block::do_add_hidden_options_to_description(options_description &arg_desc,           ///< The options_description to which the options are added
                                                                             const size_t        &/*arg_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                             ) {
	const string file_varname   { "<file>" };
	const auto sorted_links_to_file_notifier = [&] (const path &x) { the_spec.set_sorted_links_to_file( x ); };

	arg_desc.add_options()
		(
			PO_SORTED_LINKS_TO_FILE.c_str(),
			value<path>()
				->value_name   ( file_varname                  )
				->notifier     ( sorted_links_to_file_notifier ),
			(   "Rewrite the links to file " + file_varname + " (or '-' for stdout),\n"
			  + "sorted such that TCluster should cluster them identically to cath-cluster" ).c_str()
		);
}

/// \brief Generate a description of any problem that makes the specified cath_cluster_output_options_block invalid
///        or none otherwise
str_opt cath_cluster_output_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                             ) const {
	return none;
}

/// \brief Return all options names for this block
str_vec cath_cluster_output_options_block::do_get_all_options_names() const {
	return {
		cath_cluster_output_options_block::PO_CLUSTERS_TO_FILE,
		cath_cluster_output_options_block::PO_MERGES_TO_FILE,
		cath_cluster_output_options_block::PO_CLUST_SPANS_TO_FILE,
		cath_cluster_output_options_block::PO_REPS_TO_FILE,
		cath_cluster_output_options_block::PO_SORTED_LINKS_TO_FILE,
	};
}

/// \brief Getter for the cath_cluster_output_spec that the cath_cluster_output_options_block configures
const cath_cluster_output_spec & cath_cluster_output_options_block::get_cath_cluster_output_spec() const {
	return the_spec;
}
