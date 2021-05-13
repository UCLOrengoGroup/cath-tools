/// \file
/// \brief The clustmap_input_options_block class definitions

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

#include "clustmap_input_options_block.hpp"

#include <filesystem>

#include "cath/common/clone/make_uptr_clone.hpp"

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::boost::program_options::bool_switch;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::filesystem::path;
using ::std::string;
using ::std::unique_ptr;

/// \brief The option name for the cluster-membership file for the working clusters
const string clustmap_input_options_block::PO_WORKING_CLUSTMEMB_FILE  { "working-clustmemb-file"  };

/// \brief The option name for an optional file specify a cluster-membership file for map-from clusters
const string clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE { "map-from-clustmemb-file" };

/// \brief The option name for whether to read batches from working_clustmemb_file (rather than cluster membership directly)
const string clustmap_input_options_block::PO_READ_BATCHES_FROM_INPUT { "read-batches-from-input" };

/// \brief A standard do_clone method
unique_ptr<options_block> clustmap_input_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string clustmap_input_options_block::do_get_block_name() const {
	return "Input";
}

/// \brief Add this block's options to the provided options_description
void clustmap_input_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                         const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                         ) {
	const string file_varname { "<file>" };

	const auto map_from_clustmemb_file_notifier = [&] (const path_opt &x) { the_spec.set_map_from_clustmemb_file( x ); };
	const auto read_batches_from_input_notifier = [&] (const bool     &x) { the_spec.set_read_batches_from_input( x ); };

	prm_desc.add_options()
		(
			( PO_MAP_FROM_CLUSTMEMB_FILE ).c_str(),
			value<path>()
				->value_name   ( file_varname                                         )
				->notifier     ( map_from_clustmemb_file_notifier                     ),
			( "Map numbers from previous clusters specified in " + file_varname + " to their equivalents in the working clusters where possible and\n"
				"if all the cluster names in " + file_varname + " are positive integers, renumber leftover new clusters from one plus the largest\n"
				"or if not, rename with new_cmc_cluster_1, new_cmc_cluster_2, ...\n"
				"(of, if unspecified, renumber all working clusters from 1 upwards)" ).c_str()
		)
		(
			( PO_READ_BATCHES_FROM_INPUT ).c_str(),
			bool_switch()
				->notifier     ( read_batches_from_input_notifier                     )
				->default_value( clustmap_input_spec::DEFAULT_READ_BATCHES_FROM_INPUT ),
			"Read batches of work from the input file with lines of format: `batch_id working_clust_memb_file prev_clust_memb_file` where:\n"
				" * batch_id             is a unique label for the batch (with no whitespace)\n"
				" * prev_clust_memb_file is optional"
		);

	static_assert( ! clustmap_input_spec::DEFAULT_READ_BATCHES_FROM_INPUT,
		"If clustmap_input_spec::DEFAULT_READ_BATCHES_FROM_INPUT isn't false, it might mess up the bool switch in here" );
}

/// \brief Add a hidden option to the options_description for the input file
void clustmap_input_options_block::do_add_hidden_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                        const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                        ) {
	const string file_varname { "<file>" };

	const auto working_clustmemb_file_notifier = [&] (const path &x) { the_spec.set_working_clustmemb_file( x ); };

	prm_desc.add_options()
		(
			( PO_WORKING_CLUSTMEMB_FILE ).c_str(),
			value<path>()
				->value_name   ( file_varname                                         )
				->notifier     ( working_clustmemb_file_notifier                      ),
			( "Read the working cluster membership from file " + file_varname ).c_str()
		);
}

/// \brief Generate a description of any problem that makes the specified clustmap_input_options_block invalid
///        or nullopt otherwise
str_opt clustmap_input_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                        ) const {
	return get_invalid_description( get_clustmap_input_spec() );
}

/// \brief Return all options names for this block
str_vec clustmap_input_options_block::do_get_all_options_names() const {
	return {
		clustmap_input_options_block::PO_WORKING_CLUSTMEMB_FILE,
		clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE,
		clustmap_input_options_block::PO_READ_BATCHES_FROM_INPUT,
	};
}

/// \brief Getter for the clustmap_input_spec that the clustmap_input_options_block configures
const clustmap_input_spec & clustmap_input_options_block::get_clustmap_input_spec() const {
	return the_spec;
}
