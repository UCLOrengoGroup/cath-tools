/// \file
/// \brief The cath_cluster_output_spec class definitions

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

#include "cath_cluster_output_spec.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/count.hpp>

#include "common/algorithm/is_uniq.hpp"
#include "common/algorithm/sort_uniq_copy.hpp"
#include "common/cpp20/make_array.hpp"

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::std::literals::string_literals;

using ::boost::filesystem::path;
using ::boost::none;
using ::boost::numeric_cast;
using ::boost::range::count;

/// \brief Getter for an optional file to which clusters should be written
const path_opt & cath_cluster_output_spec::get_clusters_to_file() const {
	return clusters_to_file;
}

/// \brief Getter for an optional file to which merges should be written
const path_opt & cath_cluster_output_spec::get_merges_to_file() const {
	return merges_to_file;
}

/// \brief Getter for an optional file to which clust_spans should be written
const path_opt & cath_cluster_output_spec::get_clust_spans_to_file() const {
	return clust_spans_to_file;
}

/// \brief Getter for an optional file to which reps should be written
const path_opt & cath_cluster_output_spec::get_reps_to_file() const {
	return reps_to_file;
}

/// \brief Getter for an optional file to which sorted_links should be written
const path_opt & cath_cluster_output_spec::get_sorted_links_to_file() const {
	return sorted_links_to_file;
}

/// \brief Setter for an optional file to which clusters should be written
cath_cluster_output_spec & cath_cluster_output_spec::set_clusters_to_file(const path_opt &arg_clusters_to_file ///< An optional file to which clusters should be written
                                                                          ) {
	clusters_to_file = arg_clusters_to_file;
	return *this;
}

/// \brief Setter for an optional file to which merges should be written
cath_cluster_output_spec & cath_cluster_output_spec::set_merges_to_file(const path_opt &arg_merges_to_file ///< An optional file to which merges should be written
                                                                        ) {
	merges_to_file = arg_merges_to_file;
	return *this;
}

/// \brief Setter for an optional file to which clust_spans should be written
cath_cluster_output_spec & cath_cluster_output_spec::set_clust_spans_to_file(const path_opt &arg_clust_spans_to_file ///< An optional file to which clust_spans should be written
                                                                             ) {
	clust_spans_to_file = arg_clust_spans_to_file;
	return *this;
}

/// \brief Setter for an optional file to which reps should be written
cath_cluster_output_spec & cath_cluster_output_spec::set_reps_to_file(const path_opt &arg_reps_to_file ///< An optional file to which reps should be written
                                                                      ) {
	reps_to_file = arg_reps_to_file;
	return *this;
}

/// \brief Setter for an optional file to which sorted_links should be written
cath_cluster_output_spec & cath_cluster_output_spec::set_sorted_links_to_file(const path_opt &arg_sorted_links_to_file ///< An optional file to which sorted_links should be written
                                                                              ) {
	sorted_links_to_file = arg_sorted_links_to_file;
	return *this;
}

/// \brief Get the number of output paths implied by the specified cath_cluster_output_spec
///
/// \relates cath_cluster_output_spec
size_t cath::clust::get_num_output_paths(const cath_cluster_output_spec &arg_output_spec ///< The cath_cluster_output_spec to query
                                         ) {
	const auto file_presences = make_array(
		static_cast<bool>( arg_output_spec.get_clusters_to_file    () ),
		static_cast<bool>( arg_output_spec.get_merges_to_file      () ),
		static_cast<bool>( arg_output_spec.get_clust_spans_to_file () ),
		static_cast<bool>( arg_output_spec.get_reps_to_file        () ),
		static_cast<bool>( arg_output_spec.get_sorted_links_to_file() )
	);
	return static_cast<size_t>( count( file_presences, true ) );
}

/// \brief Get all output paths implied by the specified cath_cluster_output_spec
///
/// \relates cath_cluster_output_spec
path_vec cath::clust::get_all_output_paths(const cath_cluster_output_spec &arg_output_spec ///< The cath_cluster_output_spec to query
                                           ) {
	path_vec the_paths;
	if ( arg_output_spec.get_clusters_to_file() ) {
		the_paths.push_back( *arg_output_spec.get_clusters_to_file() );
	}
	if ( arg_output_spec.get_merges_to_file() ) {
		the_paths.push_back( *arg_output_spec.get_merges_to_file() );
	}
	if ( arg_output_spec.get_clust_spans_to_file() ) {
		the_paths.push_back( *arg_output_spec.get_clust_spans_to_file() );
	}
	if ( arg_output_spec.get_reps_to_file() ) {
		the_paths.push_back( *arg_output_spec.get_reps_to_file() );
	}
	if ( arg_output_spec.get_sorted_links_to_file() ) {
		the_paths.push_back( *arg_output_spec.get_sorted_links_to_file() );
	}
	return the_paths;
}

/// \brief Generate a description of any problem that makes the specified cath_cluster_output_spec invalid
///        or none otherwise
///
/// \relates cath_cluster_output_spec
str_opt cath::clust::get_invalid_description(const cath_cluster_output_spec &arg_output_spec ///< The cath_cluster_output_spec to query
                                             ) {
	const auto all_sorted_paths = sort_copy( get_all_output_paths( arg_output_spec ) );

	const size_t num_stdouts = numeric_cast<size_t>( count( all_sorted_paths, path{ "-" } ) );
	if ( num_stdouts > 1 ) {
		return "Cannot send more than one type of output to stdout (which is specified as file \"-\")"s;
	}

	if ( ! is_uniq( all_sorted_paths ) ) {
		return "Cannot send more than one type of output to the same output file"s;
	}

	return none;
}
