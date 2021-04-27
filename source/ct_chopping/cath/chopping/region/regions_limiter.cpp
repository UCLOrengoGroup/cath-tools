/// \file
/// \brief The regions_limiter class definitions

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

#include "regions_limiter.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <spdlog/spdlog.h>

#include "cath/biocore/residue_id.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/optional/make_optional_if.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::chop::detail;
using namespace ::cath::common;

using ::boost::adaptors::filtered;
using ::boost::adaptors::transformed;
using ::boost::algorithm::any_of;
using ::boost::algorithm::join;
using ::boost::make_optional;
using ::boost::none;

/// \brief Sanity check this regions_limiter and throw if there's a problem
///
/// Checks:
///  * if there are any regions without a chain label, throw
void regions_limiter::sanity_check() const {
	// If regions have been specified, check they all have chain labels
	if ( regions ) {
		if ( any_of( regions->get(), [] (const region &x) { return ! has_chain_label( x ); } ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot use regions_limiter for regions that don't have a chain_label"));
		}
	}
}

/// \brief Ctor from a vector of regions to which processing should be limited
///
/// \pre prm_regions must be non-overlapping in the context of the residue_ids they will be used to limit
regions_limiter::regions_limiter(const region_vec &prm_regions ///< The regions to which a series of residue_ids should be restricted
                                 ) : regions     { prm_regions                                            },
                                     regions_seen{ region_seen_vec( prm_regions.size(), region_seen::NO ) } {
	sanity_check();
}

/// \brief Ctor from a vector of regions to which processing should be limited
///        or none to not limit
///
/// \pre prm_regions must be non-overlapping in the context of the residue_ids they will be used to limit
regions_limiter::regions_limiter(const region_vec_opt &prm_regions
                                 ) : regions{
                                     	prm_regions
                                     	? ::boost::make_optional( std::cref( *prm_regions ) )
                                     	: none
                                     },
                                     regions_seen{
                                     	prm_regions
                                     	? ::boost::make_optional( region_seen_vec( prm_regions->size(), region_seen::NO ) )
                                     	: none
                                     } {
	sanity_check();
}

/// \brief Update state to reflect the specified residue_id as the current and return
///        whether that residue_id should be included in processing (ie if regions have
///        been specified, then this is within them)
bool regions_limiter::update_residue_is_included(const residue_id &prm_residue_id ///< The residue_id to query
                                                 ) {
	// If not restricting to any regions, then just return true
	if ( ! regions ) {
		return true;
	}

	// If currently within a region...
	if ( active_region_idx ) {
		const auto &active_region = regions->get()[ *active_region_idx ];

		// If the active region specifies start and stop
		if ( active_region.has_starts_stops() ) {

			// If prm_residue_id matches the stop of the active region, then deactivate that active region
			if ( prm_residue_id == residue_id{ get_chain_label( active_region ), get_stop_name( active_region ) } ) {
				active_region_idx = none;
			}

			// ...and either way, return true
			return true;
		}
		
		// Else the active region is a whole-chain region
		// If still inside the active whole-chain region's chain, just return true
		if ( prm_residue_id.get_chain_label() == get_chain_label( active_region ) ) {
			return true;
		}

		// Otherwise, deactivate that active region and then fall through to determine
		// whether prm_residue_id is in the start of another region
		active_region_idx = none;
	}

	// If not currently within a region, loop over the regions
	const auto num_regions = regions->get().size();
	for (const size_t &region_ctr : indices( num_regions ) ) {
		const region &the_region = regions->get()[ region_ctr ];

		// If this region specifies start and stop
		if ( the_region.has_starts_stops() ) {
			// If prm_residue_id matches the start of this region...
			if ( prm_residue_id == residue_id{ get_chain_label( the_region ), get_start_name( the_region ) } ) {
				// If prm_residue_id doesn't also match the stop of this region, then make this the active region
				if ( prm_residue_id != residue_id{ get_chain_label( the_region ), get_stop_name( the_region ) } ) {
					active_region_idx = region_ctr;
				}

				// ...and either way, mark the region as seen and return true
				regions_seen.get()[ region_ctr ] = region_seen::YES;
				return true;
			}
		}
		// Else this region is a whole-chain region
		else {
			// If inside this whole-chain region's chain, then make this the active region,
			// mark the region as seen and return true
			if ( prm_residue_id.get_chain_label() == get_chain_label( the_region ) ) {
				active_region_idx = region_ctr;
				regions_seen.get()[ region_ctr ] = region_seen::YES;
				return true;
			}
		}
	}

	// Otherwise, not inside any region so return false
	return false;
}

/// \brief Return any specified regions that remain unseen
region_vec regions_limiter::unseen_regions() const {
	// If not restricting to any regions, then just return the empty results
	if ( ! regions ) {
		return {};
	}

	// Otherwise, build a vector of any regions that haven't been seen
	return transform_build<region_vec>(
		indices( regions->get().size() )
			| filtered( [&] (const size_t &region_idx) {
				return ( regions_seen.get()[ region_idx ] == region_seen::NO );
			} ),
		[&] (const size_t &region_idx) {
			return regions->get()[ region_idx ];
		}
	);
}

/// \brief Return a string describing any of the specified regions_limiter's specified regions that it hasn't yet seen
///        or return none if no regions were specified or all have been seen
///
/// \relates regions_limiter
str_opt cath::chop::warn_str_if_specified_regions_remain_unseen(const regions_limiter &prm_regions ///< The regions_limiter to query
                                                                ) {
	const auto unseen_regions = prm_regions.unseen_regions();
	return make_optional_if_fn(
		! unseen_regions.empty(),
		[&] {
			return
				"Unable to find anything corresponding to region(s) : "
				+ join(
					unseen_regions
						| transformed( [] (const region &x) { return to_string( x ); } ),
					", "
				);
		}
	);
}

/// \brief Warn about any of the specified regions_limiter's specified regions that it hasn't yet seen
///        or do nothing if no regions were specified or all have been seen
///
/// \relates regions_limiter
void cath::chop::warn_if_specified_regions_remain_unseen(const regions_limiter &prm_regions ///< The regions_limiter to query
                                                         ) {
	const auto warn_str = warn_str_if_specified_regions_remain_unseen( prm_regions );
	if ( warn_str ) {
		::spdlog::warn( *warn_str );
	}
}
