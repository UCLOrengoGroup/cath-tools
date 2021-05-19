/// \file
/// \brief The cath_cluster_clustering_options_block class definitions

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

#include "cath_cluster_clustering_options_block.hpp"

#include "cath/cath_cluster/options/spec/clustering_levels.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::nullopt;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<options_block> cath_cluster_clustering_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string cath_cluster_clustering_options_block::do_get_block_name() const {
	return "Clustering";
}

/// \brief Add this block's options to the provided options_description
void cath_cluster_clustering_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                                  const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                                  ) {
	const string levels_varname { "<levels>" };

	const auto levels_notifier = [&] (const clustering_levels &x) { the_spec.set_levels( x.levels ); };

	prm_desc.add_options()
		(
			string( PO_LEVELS ).c_str(),
			value<clustering_levels>()
				->value_name   ( levels_varname  )
				->notifier     ( levels_notifier )
				->required     (                 ),
			( "Cluster at levels "
			  + levels_varname
			  + ", which is ordered values separated by commas (eg 35,60,95,100)" ).c_str()
		);
}

/// \brief Generate a description of any problem that makes the specified cath_cluster_clustering_options_block invalid
///        or nullopt otherwise
str_opt cath_cluster_clustering_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                                 ) const {
	return nullopt;
}

/// \brief Return all options names for this block
str_view_vec cath_cluster_clustering_options_block::do_get_all_options_names() const {
	return {
		PO_LEVELS,
	};
}

/// \brief Getter for the cath_cluster_clustering_spec that the cath_cluster_clustering_options_block configures
const cath_cluster_clustering_spec & cath_cluster_clustering_options_block::get_cath_cluster_clustering_spec() const {
	return the_spec;
}
