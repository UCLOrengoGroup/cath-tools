/// \file
/// \brief The full_hit_list class definitions

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

#include "full_hit_list.h"

#include <boost/range/combine.hpp>

#include "common/algorithm/append.h"
#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/range/max_proj_element.h"
#include "resolve_hits/trim/hit_seg_boundary_fns.h"

using namespace cath::common;
using namespace cath::rslv;

using boost::make_optional;
using boost::none;
using boost::range::combine;
using std::ostream;
using std::string;

/// \brief Get a vector of all the segments in the specified list of hits that
///        differ from the specified hit
///
/// \relates full_hit_list
hit_seg_vec get_other_hits_segments(const full_hit      &arg_full_hit, ///< The hit whose segments should be excluded from the list
                                    const full_hit_list &arg_full_hits ///< The list of hits from which the segments should be drawn
                                    ) {
	hit_seg_vec results;
	for (const full_hit &the_full_hit : arg_full_hits) {
		if ( the_full_hit != arg_full_hit ) {
			append( results, the_full_hit.get_segments() );
		}
	}
	return results;
}

/// \brief Calculate the resolved boundaries for the specified hit in the context of
///        the specified list of hits and trim_spec
///
/// \relates full_hit_list
seg_boundary_pair_vec cath::rslv::resolved_boundaries(const full_hit      &arg_full_hit,  ///< The hit whose boundaries should be resolved
                                                      const full_hit_list &arg_full_hits, ///< The list of hits providing the context (may include the original hit)
                                                      const trim_spec     &arg_trim_spec  ///< The trim_spec to use
                                                      ) {
	const hit_seg_vec other_segments = get_other_hits_segments( arg_full_hit, arg_full_hits );

	// Build up a vector of res_arr_res_arr_pairs, one calculated for each segment of arg_full_hit
	return transform_build<seg_boundary_pair_vec>(
		arg_full_hit.get_segments(),
		[&] (const hit_seg &x) {
			return get_boundary_pair( x, other_segments, arg_trim_spec );
		}
	);
}

/// \brief Merge original boundaries with resolved boundaries, replacing originals with
///        resolveds where they exist
hit_seg_vec cath::rslv::merge_boundaries(const hit_seg_vec           &arg_segs,               ///< The original boundaries
                                         const seg_boundary_pair_vec &arg_resolved_boundaries ///< The resolved boundaries (where resolving has been required)
                                         ) {
	return transform_build<hit_seg_vec>(
		combine( arg_segs, arg_resolved_boundaries ),
		[] (const boost::tuple<hit_seg, seg_boundary_pair> &x) {
			return hit_seg{
				x.get<1>().first .value_or( x.get<0>().get_start_arrow() ),
				x.get<1>().second.value_or( x.get<0>().get_stop_arrow () )
			};
		}
	);
}

/// \brief Resolved all boundaries for the specified hit in the context of
///        the specified list of hits and trim_spec (ie return either the
///        resolved boundary or the original where no resolving is required)
///
/// \relates full_hit_list
hit_seg_vec cath::rslv::resolve_all_boundaries(const full_hit      &arg_full_hit,  ///< The hit whose boundaries should be resolved
                                               const full_hit_list &arg_full_hits, ///< The list of hits providing the context (may include the original hit)
                                               const trim_spec     &arg_trim_spec  ///< The trim_spec to use
                                               ) {
	return merge_boundaries(
		arg_full_hit.get_segments(),
		resolved_boundaries( arg_full_hit, arg_full_hits, arg_trim_spec )
	);
}

/// \brief Get a string describing the full resolved boundaries for the specified hit in the context of
///        the specified list of hits and trim_spec
///
/// \relates full_hit_list
string cath::rslv::get_all_resolved_segments_string(const full_hit      &arg_full_hit,  ///< The hit whose boundaries should be resolved
                                                    const full_hit_list &arg_full_hits, ///< The list of hits providing the context (may include the original hit)
                                                    const trim_spec     &arg_trim_spec  ///< The trim_spec to use
                                                    ) {
	return get_segments_string(
		resolve_all_boundaries( arg_full_hit, arg_full_hits, arg_trim_spec ),
		none
	);
}

/// \brief Generate a string describing the specified full_hit_list
///
/// \relates full_hit_list
string cath::rslv::to_string(const full_hit_list &arg_full_hit_list ///< The full_hit_list to describe
                             ) {
	return "full_hit_list["
		+ ::std::to_string( arg_full_hit_list.size() )
		+ "full_hits]";
}

/// \brief Insert a description of the specified full_hit_list into the specified ostream
///
/// \relates full_hit_list
ostream & cath::rslv::operator<<(ostream             &arg_ostream, ///< The ostream into which the description should be inserted
                                 const full_hit_list &arg_full_hit_list ///< The full_hit_list to describe
                                 ) {
	arg_ostream << to_string( arg_full_hit_list );
	return arg_ostream;
}

/// \brief Get the maximum stop residue of all the full_hits in the specified full_hit_list
///        or none if the full_hit_list is empty
///
/// \relates full_hit_list
residx_opt cath::rslv::get_max_stop(const full_hit_list &arg_full_hit_list ///< The full_hit_list to query
                                    ) {
	// Can't just use make_optional(bool, T) because the second part shouldn't be evaluated if arg_full_hit_list's empty
	return arg_full_hit_list.empty()
		? none
		: make_optional( max_proj(
			arg_full_hit_list,
			std::less<>{},
			&cath::rslv::get_stop_res_arrow
		).res_before() );
}


/// \brief Get the best (ie maximum) score of all the full_hits in the specified full_hit_list
///        or none if the full_hit_list is empty
///
/// \relates full_hit_list
resscr_opt cath::rslv::get_best_crh_score(const full_hit_list  &arg_full_hit_list, ///< The full_hit_list to query
                                          const crh_score_spec &arg_score_spec     ///< The score_spec specifying how to calculate the crh-score of a full-hit
                                          ) {
	// Can't just use make_optional(bool, T) because the second part shouldn't be evaluated if arg_full_hit_list's empty
	return arg_full_hit_list.empty()
		? none
		: make_optional( max_proj(
			arg_full_hit_list,
			std::less<>{},
			[&] (const full_hit &x) { return get_crh_score( x, arg_score_spec ); }
		) );
}
