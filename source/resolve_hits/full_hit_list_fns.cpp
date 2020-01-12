/// \file
/// \brief The full_hit_list functions definitions

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

#include "full_hit_list_fns.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/combine.hpp>

#include "common/algorithm/append.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/max_proj_element.hpp"
#include "resolve_hits/file/alnd_rgn.hpp"
#include "resolve_hits/full_hit_fns.hpp"
#include "resolve_hits/full_hit_rapidjson.hpp"
#include "resolve_hits/trim/seq_seg_boundary_fns.hpp"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::seq;

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::make_optional;
using boost::none;
using boost::range::combine;
using std::string;


/// \brief Generate a string describing the specified full_hit_list in the specified format
///
/// This is separated from a normal to_string() function because
/// the interface requirements of this may change (eg to demand that the client
/// passes the crh_score_spec, crh_segment_spec, hits_boundary_output)
///
/// \relates full_hit_list
string cath::rslv::to_output_string(const full_hit_list       &prm_full_hits,        ///< The list of full_hits to describe
                                    const crh_segment_spec    &prm_crh_segment_spec, ///< The crh_segment_spec specifying any trimming that should be performed on the output segments
                                    const hit_output_format   &prm_format,           ///< The format in which the hit_arch should be described
                                    const string              &prm_prefix,           ///< Any prefix that should come before the hit in hit_output_format::JON
                                    const hit_boundary_output &prm_boundary_output   ///< Whether to output the trimmed or original boundaries
                                    ) {
	const bool is_jon = ( prm_format == hit_output_format::JON );
	const string prefix    = is_jon ? ""   : "full_hit_list[\n\t";
	const string separator = is_jon ? "\n" : "\n\t";
	const string suffix    = is_jon ? "\n" : "\n]";

	return prefix
		+ join(
			prm_full_hits
				| transformed( [&] (const full_hit &x) {
					return to_string(
						x,
						prm_format,
						prm_prefix,
						make_optional( prm_crh_segment_spec ),
						make_optional( prm_full_hits ),
						prm_boundary_output
					);
				} ),
			separator
		)
		+ suffix;
}

/// \brief Get the best (ie maximum) score of all the full_hits in the specified full_hit_list
///        or none if the full_hit_list is empty
///
/// \relates full_hit_list
resscr_opt cath::rslv::get_best_crh_score(const full_hit_list  &prm_full_hit_list, ///< The full_hit_list to query
                                          const crh_score_spec &prm_score_spec     ///< The score_spec specifying how to calculate the crh-score of a full-hit
                                          ) {
	// Can't just use make_optional(bool, T) because the second part shouldn't be evaluated if prm_full_hit_list's empty
	return prm_full_hit_list.empty()
		? none
		: make_optional( max_proj(
			prm_full_hit_list,
			std::less<>{},
			[&] (const full_hit &x) { return get_crh_score( x, prm_score_spec ); }
		) );
}

/// \brief Get a vector of all the segments in the specified list of hits that
///        differ from the specified hit
///
/// \relates full_hit_list
static seq_seg_vec get_other_hits_segments(const full_hit      &prm_full_hit, ///< The hit whose segments should be excluded from the list
                                           const full_hit_list &prm_full_hits ///< The list of hits from which the segments should be drawn
                                           ) {
	seq_seg_vec results;
	for (const full_hit &the_full_hit : prm_full_hits) {
		if ( the_full_hit != prm_full_hit ) {
			append( results, the_full_hit.get_segments() );
		}
	}
	return results;
}

/// \brief Calculate the resolved boundaries for the specified hit in the context of
///        the specified list of hits and crh_segment_spec
///
/// \relates full_hit_list
seg_boundary_pair_vec cath::rslv::resolved_boundaries(const full_hit         &prm_full_hit,        ///< The hit whose boundaries should be resolved
                                                      const full_hit_list    &prm_full_hits,       ///< The list of hits providing the context (may include the original hit)
                                                      const crh_segment_spec &prm_crh_segment_spec ///< The crh_segment_spec to use
                                                      ) {
	const seq_seg_vec other_segments = get_other_hits_segments( prm_full_hit, prm_full_hits );

	// Build up a vector of res_arr_res_arr_pairs, one calculated for each segment of prm_full_hit
	return transform_build<seg_boundary_pair_vec>(
		prm_full_hit.get_segments(),
		[&] (const seq_seg &x) {
			return get_boundary_pair( x, other_segments, prm_crh_segment_spec );
		}
	);
}

/// \brief Merge original boundaries with resolved boundaries, replacing originals with
///        resolveds where they exist
seq_seg_opt_vec cath::rslv::merge_boundaries(const seq_seg_vec           &prm_segs,                ///< The original boundaries
                                             const seg_boundary_pair_vec &prm_resolved_boundaries, ///< The resolved boundaries (where resolving has been required)
                                             const crh_segment_spec      &prm_crh_segment_spec     ///< The crh_segment_spec defining how the segments will be handled (eg trimmed) by the algorithm
                                             ) {
	return transform_build<seq_seg_opt_vec>(
		combine( prm_segs, prm_resolved_boundaries ),
		[&] (const boost::tuple<seq_seg, seg_boundary_pair> &x) {
			const auto &the_seq_seg = x.get<0>();
			const auto &seg_bnd_1   = x.get<1>().first;
			const auto &seg_bnd_2   = x.get<1>().second;
			return make_optional(
				get_length( the_seq_seg ) >= prm_crh_segment_spec.get_min_seg_length(),
				seq_seg{
					seg_bnd_1.value_or( the_seq_seg.get_start_arrow() ),
					seg_bnd_2.value_or( the_seq_seg.get_stop_arrow () )
				}
			);
		}
	);
}

/// \brief Resolved all boundaries for the specified hit in the context of
///        the specified list of hits and crh_segment_spec (ie return either the
///        resolved boundary or the original where no resolving is required)
///
/// \relates full_hit_list
seq_seg_opt_vec cath::rslv::resolve_all_boundaries(const full_hit         &prm_full_hit,        ///< The hit whose boundaries should be resolved
                                                   const full_hit_list    &prm_full_hits,       ///< The list of hits providing the context (may include the original hit)
                                                   const crh_segment_spec &prm_crh_segment_spec ///< The crh_segment_spec to use
                                                   ) {
	return merge_boundaries(
		prm_full_hit.get_segments(),
		resolved_boundaries( prm_full_hit, prm_full_hits, prm_crh_segment_spec ),
		prm_crh_segment_spec
	);
}

/// \brief Get a string describing the full resolved boundaries for the specified hit in the context of
///        the specified list of hits and crh_segment_spec
///
/// \relates full_hit_list
string cath::rslv::get_all_resolved_segments_string(const full_hit         &prm_full_hit,        ///< The hit whose boundaries should be resolved
                                                    const full_hit_list    &prm_full_hits,       ///< The list of hits providing the context (may include the original hit)
                                                    const crh_segment_spec &prm_crh_segment_spec ///< The crh_segment_spec to use
                                                    ) {
	return get_segments_string(
		resolve_all_boundaries( prm_full_hit, prm_full_hits, prm_crh_segment_spec ),
		none
	);
}

/// \brief Write the specified full_hit_list to the specified rapidjson_writer
string cath::rslv::to_json_string_with_compact_fullhits(const full_hit_list        &prm_full_hit_list, ///< The full_hit_list to write
                                                        const crh_segment_spec_opt &prm_segment_spec,  ///< An optional crh_segment_spec which can be used for including each full_hit's trimmed boundaries and resolved boundaries
                                                        const size_t               &prm_extra_depth    ///< The number of levels of depth
                                                        ) {
	return string_of_rapidjson_write<json_style::PRETTY> (
		[&] (rapidjson_writer<json_style::PRETTY> &x) {
			write_to_rapidjson_with_compact_fullhits( x, prm_full_hit_list, prm_segment_spec );
		},
		prm_extra_depth
	);
}

