/// \file
/// \brief The cath_cluster_clustering_spec class definitions

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

#include "cath_cluster_clustering_spec.hpp"

#include <boost/algorithm/cxx11/is_sorted.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/sort.hpp>

#include "common/boost_addenda/range/max_proj_element.hpp"
#include "common/exception/invalid_argument_exception.hpp"

#include <functional>

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::std::literals::string_literals;

using ::boost::adaptors::transformed;
using ::boost::algorithm::is_sorted;
using ::boost::algorithm::join;
using ::boost::algorithm::to_lower_copy;
using ::boost::lexical_cast;
using ::boost::make_optional;
using ::boost::none;
using ::boost::range::sort;
using ::std::greater;
using ::std::less;
using ::std::string;

/// \brief Getter for the levels at which the clustering should be performed
const strength_vec & cath_cluster_clustering_spec::get_levels() const {
	return levels;
}

/// \brief Setter for the levels at which the clustering should be performed
cath_cluster_clustering_spec & cath_cluster_clustering_spec::set_levels(const strength_vec &prm_levels ///< The levels at which the clustering should be performed
                                                                                  ) {
	levels = prm_levels;
	return *this;
}

/// \brief Generate a description of any problem that makes the specified cath_cluster_clustering_spec invalid
///        or none otherwise
///
/// \relates cath_cluster_clustering_spec
str_opt cath::clust::get_invalid_description(const cath_cluster_clustering_spec &prm_clustering_spec ///< The cath_cluster_clustering_spec to query
                                             ) {
	if ( prm_clustering_spec.get_levels().empty() ) {
		return "Most specify at least one clustering level"s;
	}

	return none;
}

/// \brief Get a warning string if the specified clustering levels aren't sorted suitable for the specified link_dirn
///        or none otherwise
str_opt cath::clust::get_dissim_sort_warning(const strength_vec &prm_levels,   ///< The clustering levels to check
                                             const link_dirn    &prm_link_dirn ///< Whether the links in the input file represent strengths or dissimilarities
                                             ) {
	const bool is_correct = ( prm_link_dirn == link_dirn::STRENGTH )
		? is_sorted( prm_levels, less   <>{} )
		: is_sorted( prm_levels, greater<>{} );
	return make_optional(
		! is_correct,
		"The levels ("
			+ join(
				prm_levels
					| transformed( [] (const strength &x) { return lexical_cast<string>( x ); } ),
				", "
			)
			+ ") are not sorted to be "
			+ ( ( prm_link_dirn == link_dirn::STRENGTH ) ? "increasing" : "decreasing" )
			+ " as would be expected with a "
			+ to_lower_copy( to_string( prm_link_dirn ) )
			+ " link direction"
	);
}

/// \brief Convert the specified clustering levels into dissimilarities
///        (as necessary according to the specified link_dirn) and sort
void cath::clust::make_dissim_and_sort(strength_vec    &prm_levels,   ///< The clustering levels to alter
                                       const link_dirn &prm_link_dirn ///< Whether the links in the input file represent strengths or dissimilarities
                                       ) {
	if ( prm_link_dirn == link_dirn::STRENGTH ) {
		for (strength &x : prm_levels) {
			x = -x;
		}
	}
	sort( prm_levels );
}

/// \brief Convert a copy of the specified clustering levels into dissimilarities
///        (as necessary according to the specified link_dirn), sort and return
strength_vec cath::clust::make_dissim_and_sort_copy(strength_vec     prm_levels,   ///< The clustering levels to query
                                                    const link_dirn &prm_link_dirn ///< Whether the links in the input file represent strengths or dissimilarities
                                                    ) {
	make_dissim_and_sort( prm_levels, prm_link_dirn );
	return prm_levels;
}

/// \brief Get the maximum dissimilarity value corresponding to the specified clustering levels which are in the specified link_dirn
strength cath::clust::get_max_dissim(const strength_vec &prm_levels,   ///< The clustering levels to query
                                     const link_dirn    &prm_link_dirn ///< Whether the links in the input file represent strengths or dissimilarities
                                     ) {
	if ( prm_levels.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get the max dissimilarity of an empty list of levels"));
	}
	return max_proj(
		prm_levels,
		std::less<>{},
		[&] (const strength &x) { return ( prm_link_dirn == link_dirn::STRENGTH ) ? -x : x; }
	);
}

/// \brief Get a warning string if the specified cath_cluster_clustering_spec's clustering levels aren't
///        sorted suitable for the specified link_dirn or none otherwise
///
/// \relates cath_cluster_clustering_spec
str_opt cath::clust::get_dissim_sort_warning(const cath_cluster_clustering_spec &prm_spec,     ///< The cath_cluster_clustering_spec to check
                                             const link_dirn                    &prm_link_dirn ///< Whether the links in the input file represent strengths or dissimilarities
                                             ) {
	return get_dissim_sort_warning( prm_spec.get_levels(), prm_link_dirn );
}

/// \brief Convert a copy of the specified cath_cluster_clustering_spec's clustering levels into dissimilarities
///        (as necessary according to the specified link_dirn), sort and return
///
/// \relates cath_cluster_clustering_spec
strength_vec cath::clust::get_sorted_dissims(const cath_cluster_clustering_spec &prm_spec,     ///< The cath_cluster_clustering_spec to query
                                             const link_dirn                    &prm_link_dirn ///< Whether the links in the input file represent strengths or dissimilarities
                                             ) {
	return make_dissim_and_sort_copy( prm_spec.get_levels(), prm_link_dirn );
}

/// \brief Get the maximum dissimilarity value corresponding to the specified cath_cluster_clustering_spec's
///        clustering levels which are in the specified link_dirn
///
/// \relates cath_cluster_clustering_spec
strength cath::clust::get_max_dissim(const cath_cluster_clustering_spec &prm_spec,     ///< The cath_cluster_clustering_spec to query
                                     const link_dirn                    &prm_link_dirn ///< Whether the links in the input file represent strengths or dissimilarities
                                     ) {
	return get_max_dissim( prm_spec.get_levels(), prm_link_dirn );
}
