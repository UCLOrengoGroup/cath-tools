/// \file
/// \brief The hit_seg_boundary_fns class definitions

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

#include "hit_seg_boundary_fns.hpp"

#include <boost/algorithm/cxx11/none_of.hpp>
#include <boost/optional.hpp>
#include <boost/range/adaptor/filtered.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/boost_addenda/range/max_proj_element.hpp"
#include "resolve_hits/hit_seg.hpp"
#include "resolve_hits/options/spec/crh_segment_spec.hpp"
#include "resolve_hits/res_arrow.hpp"
#include "resolve_hits/trim/resolve_boundary.hpp"
#include "resolve_hits/trim/trim_spec.hpp"

using namespace cath::common;
using namespace cath::rslv;

using boost::adaptors::filtered;
using boost::algorithm::none_of;
using boost::none;
using std::make_pair;

/// \brief Get the specified boundary for the specified segment in the context of the other segments and crh_segment_spec
///        or none if there are no overlapping segments on the relevant side of the query segment
///
/// The other segments may include the original segment (though that may
/// suggest that the client code isn't correctly handling hits)
res_arrow_opt cath::rslv::detail::get_boundary_impl(const boundary_wanted  &arg_dirn,            ///< The direction of the boundary to query
                                                    const hit_seg          &arg_segment,         ///< The hit to query
                                                    const hit_seg_vec      &arg_other_segments,  ///< The other segments from other hits
                                                    const crh_segment_spec &arg_crh_segment_spec ///< The crh_segment_spec to use
                                                    ) {
	// Make a function that returns whether the specified segment is long enough
	const auto is_long_enough_fn = [&] (const hit_seg &x) {
		return get_length( x ) >= arg_crh_segment_spec.get_min_seg_length();
	};

	// If the query segment isn't long enough then return none
	if ( ! is_long_enough_fn( arg_segment ) ) {
		return none;
	}

	// Create a range that excludes the query segment
	const auto other_segments = arg_other_segments
		| filtered( [&] (const hit_seg &x) { return x != arg_segment; } );

	// Make a function that returns whether the first seg is more in the
	// direction of the end of interest than the second
	const auto in_wanted_dirn_fn = [&] (const hit_seg &x, const hit_seg &y) {
		return ( arg_dirn == boundary_wanted::START )
			? ( midpoint_less( x, y ) )
			: ( midpoint_less( y, x ) );
	};
	// Make a function that returns whether the specified segment is
	// long enough and on the correct side of the query segment
	const auto is_long_enough_and_on_correct_side_fn = [&] (const hit_seg &x) {
		return (
			in_wanted_dirn_fn( x, arg_segment )
			&&
			is_long_enough_fn( x )
		);
	};

	// If none of the other segments are on the correct side of the query, return none
	if ( none_of( other_segments, is_long_enough_and_on_correct_side_fn ) ) {
		return none;
	}

	// Build a vector of the segments that are on the correct side and are long enough
	//
	// \todo Come a better filtered that doesn't require its fn be copyable, don't
	//       bother doing this step to build a vector, just use the filtered
	//       other_segments directly in the call to max_proj()
	const auto other_segs_on_correct_side = copy_build<hit_seg_vec>(
		other_segments | filtered( is_long_enough_and_on_correct_side_fn )
	);

	// Calculate and return the resolved boundary between the closest segment and the query
	return calc_resolved_boundary(
		max_proj(
			other_segs_on_correct_side,
			in_wanted_dirn_fn
		),
		arg_segment,
		arg_crh_segment_spec.get_overlap_trim_spec()
	);
}

/// \brief Get the start boundary for the specified segment in the context of the other segments and crh_segment_spec
res_arrow_opt cath::rslv::get_start_boundary(const hit_seg          &arg_segment,         ///< The hit to query
                                             const hit_seg_vec      &arg_other_segments,  ///< The other segments from other hits
                                             const crh_segment_spec &arg_crh_segment_spec ///< The crh_segment_spec to use
                                             ) {
	return get_boundary_impl(
		detail::boundary_wanted::START,
		arg_segment,
		arg_other_segments,
		arg_crh_segment_spec
	);
}

/// \brief Get the stop boundary for the specified segment in the context of the other segments and crh_segment_spec
res_arrow_opt cath::rslv::get_stop_boundary(const hit_seg          &arg_segment,         ///< The hit to query
                                            const hit_seg_vec      &arg_other_segments,  ///< The other segments from other hits
                                            const crh_segment_spec &arg_crh_segment_spec ///< The crh_segment_spec to use
                                            ) {
	return get_boundary_impl(
		detail::boundary_wanted::STOP,
		arg_segment,
		arg_other_segments,
		arg_crh_segment_spec
	);
}

/// \brief Get the start/stop boundary pair for the specified segment in the context of the other segments and crh_segment_spec
seg_boundary_pair cath::rslv::get_boundary_pair(const hit_seg          &arg_segment,        ///< The hit to query
                                                const hit_seg_vec      &arg_other_segments, ///< The other segments from other hits
                                                const crh_segment_spec &arg_crh_segment_spec       ///< The crh_segment_spec to use
                                                ) {
	return make_pair(
		get_start_boundary( arg_segment, arg_other_segments, arg_crh_segment_spec ),
		get_stop_boundary ( arg_segment, arg_other_segments, arg_crh_segment_spec )
	);
}

/// \brief Calculate shared boundary between the two specified segments given the specified trim_spec
///
/// If the segments don't overlap, then this returns none
///
/// The segments do not need to be provided in order
res_arrow_opt cath::rslv::calc_resolved_boundary(const hit_seg   &arg_segment_a, ///< The first  segment to compare
                                                 const hit_seg   &arg_segment_b, ///< The second segment to compare
                                                 const trim_spec &arg_trim_spec  ///< The trim_spec to use
                                                 ) {
	// Get the correct order for the segments
	const bool a_is_first = midpoint_less( arg_segment_a, arg_segment_b );
	const hit_seg &seg_lhs =     a_is_first   ? arg_segment_a : arg_segment_b;
	const hit_seg &seg_rhs = ( ! a_is_first ) ? arg_segment_a : arg_segment_b;

	// If the segments don't overlap, return none
	if ( seg_lhs.get_stop_arrow() <= seg_rhs.get_start_arrow() ) {
		return none;
	}

	// Otherwise, grab the relevant trim lengths and resolve the boundary
	const auto lhs_stop_trim  = stop_trimming_of_length ( arg_trim_spec, get_length( seg_lhs ) );
	const auto rhs_start_trim = start_trimming_of_length( arg_trim_spec, get_length( seg_rhs ) );
	return resolve_boundary(
		seg_lhs.get_stop_arrow () - lhs_stop_trim,
		lhs_stop_trim,
		seg_rhs.get_start_arrow() + rhs_start_trim,
		rhs_start_trim
	);
}
