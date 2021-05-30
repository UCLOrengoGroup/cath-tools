/// \file
/// \brief The seq_seg_boundary_fns class definitions

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

#include "seq_seg_boundary_fns.hpp"

#include <boost/algorithm/cxx11/none_of.hpp>
#include <boost/range/adaptor/filtered.hpp>

#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/boost_addenda/range/max_proj_element.hpp"
#include "cath/resolve_hits/options/spec/crh_segment_spec.hpp"
#include "cath/resolve_hits/trim/resolve_boundary.hpp"
#include "cath/resolve_hits/trim/trim_spec.hpp"
#include "cath/seq/seq_arrow.hpp"
#include "cath/seq/seq_seg.hpp"

using namespace ::cath::common;
using namespace ::cath::rslv;
using namespace ::cath::seq;

using ::boost::adaptors::filtered;
using ::boost::algorithm::none_of;
using ::std::make_pair;
using ::std::nullopt;

/// \brief Get the specified boundary for the specified segment in the context of the other segments and crh_segment_spec
///        or nullopt if there are no overlapping segments on the relevant side of the query segment
///
/// The other segments may include the original segment (though that may
/// suggest that the client code isn't correctly handling hits)
res_arrow_opt cath::rslv::detail::get_boundary_impl(const boundary_wanted  &prm_dirn,            ///< The direction of the boundary to query
                                                    const seq_seg          &prm_segment,         ///< The hit to query
                                                    const seq_seg_vec      &prm_other_segments,  ///< The other segments from other hits
                                                    const crh_segment_spec &prm_crh_segment_spec ///< The crh_segment_spec to use
                                                    ) {
	// Make a function that returns whether the specified segment is long enough
	const auto is_long_enough_fn = [&] (const seq_seg &x) {
		return get_length( x ) >= prm_crh_segment_spec.get_min_seg_length();
	};

	// If the query segment isn't long enough then return nullopt
	if ( ! is_long_enough_fn( prm_segment ) ) {
		return nullopt;
	}

	// Create a range that excludes the query segment
	const auto other_segments = prm_other_segments
		| filtered( [&] (const seq_seg &x) { return x != prm_segment; } );

	// Make a function that returns whether the first seg is more in the
	// direction of the end of interest than the second
	const auto in_wanted_dirn_fn = [&] (const seq_seg &x, const seq_seg &y) {
		return ( prm_dirn == boundary_wanted::START )
			? ( midpoint_less( x, y ) )
			: ( midpoint_less( y, x ) );
	};
	// Make a function that returns whether the specified segment is
	// long enough and on the correct side of the query segment
	const auto is_long_enough_and_on_correct_side_fn = [&] (const seq_seg &x) {
		return (
			in_wanted_dirn_fn( x, prm_segment )
			&&
			is_long_enough_fn( x )
		);
	};

	// If none of the other segments are on the correct side of the query, return nullopt
	if ( none_of( other_segments, is_long_enough_and_on_correct_side_fn ) ) {
		return nullopt;
	}

	// Calculate and return the resolved boundary between the closest segment and the query
	return calc_resolved_boundary(
		max_proj(
			other_segments | filtered( is_long_enough_and_on_correct_side_fn ),
			in_wanted_dirn_fn
		),
		prm_segment,
		prm_crh_segment_spec.get_overlap_trim_spec()
	);
}

/// \brief Get the start boundary for the specified segment in the context of the other segments and crh_segment_spec
res_arrow_opt cath::rslv::get_start_boundary(const seq_seg          &prm_segment,         ///< The hit to query
                                             const seq_seg_vec      &prm_other_segments,  ///< The other segments from other hits
                                             const crh_segment_spec &prm_crh_segment_spec ///< The crh_segment_spec to use
                                             ) {
	return get_boundary_impl(
		detail::boundary_wanted::START,
		prm_segment,
		prm_other_segments,
		prm_crh_segment_spec
	);
}

/// \brief Get the stop boundary for the specified segment in the context of the other segments and crh_segment_spec
res_arrow_opt cath::rslv::get_stop_boundary(const seq_seg          &prm_segment,         ///< The hit to query
                                            const seq_seg_vec      &prm_other_segments,  ///< The other segments from other hits
                                            const crh_segment_spec &prm_crh_segment_spec ///< The crh_segment_spec to use
                                            ) {
	return get_boundary_impl(
		detail::boundary_wanted::STOP,
		prm_segment,
		prm_other_segments,
		prm_crh_segment_spec
	);
}

/// \brief Get the start/stop boundary pair for the specified segment in the context of the other segments and crh_segment_spec
seg_boundary_pair cath::rslv::get_boundary_pair(const seq_seg          &prm_segment,        ///< The hit to query
                                                const seq_seg_vec      &prm_other_segments, ///< The other segments from other hits
                                                const crh_segment_spec &prm_crh_segment_spec       ///< The crh_segment_spec to use
                                                ) {
	return make_pair(
		get_start_boundary( prm_segment, prm_other_segments, prm_crh_segment_spec ),
		get_stop_boundary ( prm_segment, prm_other_segments, prm_crh_segment_spec )
	);
}

/// \brief Calculate shared boundary between the two specified segments given the specified trim_spec
///
/// If the segments don't overlap, then this returns nullopt
///
/// The segments do not need to be provided in order
res_arrow_opt cath::rslv::calc_resolved_boundary(const seq_seg   &prm_segment_a, ///< The first  segment to compare
                                                 const seq_seg   &prm_segment_b, ///< The second segment to compare
                                                 const trim_spec &prm_trim_spec  ///< The trim_spec to use
                                                 ) {
	// Get the correct order for the segments
	const bool a_is_first = midpoint_less( prm_segment_a, prm_segment_b );
	const seq_seg &seg_lhs =     a_is_first   ? prm_segment_a : prm_segment_b;
	const seq_seg &seg_rhs = ( ! a_is_first ) ? prm_segment_a : prm_segment_b;

	// If the segments don't overlap, return nullopt
	if ( seg_lhs.get_stop_arrow() <= seg_rhs.get_start_arrow() ) {
		return nullopt;
	}

	// Otherwise, grab the relevant trim lengths and resolve the boundary
	const auto lhs_stop_trim  = stop_trimming_of_length ( prm_trim_spec, get_length( seg_lhs ) );
	const auto rhs_start_trim = start_trimming_of_length( prm_trim_spec, get_length( seg_rhs ) );
	return resolve_boundary(
		seg_lhs.get_stop_arrow () - lhs_stop_trim,
		lhs_stop_trim,
		seg_rhs.get_start_arrow() + rhs_start_trim,
		rhs_start_trim
	);
}
