/// \file
/// \brief The full_hit_list functions header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FULL_HIT_LIST_FNS_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FULL_HIT_LIST_FNS_HPP

#include <optional>

#include "cath/resolve_hits/full_hit_list.hpp"
#include "cath/resolve_hits/hit_output_format.hpp"
#include "cath/resolve_hits/options/spec/crh_score_spec.hpp"
#include "cath/resolve_hits/options/spec/crh_segment_spec.hpp"
#include "cath/resolve_hits/options/spec/hit_boundary_output.hpp"

// clang-format off
namespace cath::rslv { class read_and_process_mgr; }
// clang-format on

namespace cath::rslv {

	std::string to_output_string(const full_hit_list &,
	                             const crh_segment_spec &,
	                             const hit_output_format & = hit_output_format::CLASS,
	                             const std::string & = std::string{},
	                             const hit_boundary_output & = hit_boundary_output::ORIG);
	resscr_opt get_best_crh_score(const full_hit_list &,
	                              const crh_score_spec &);

	seg_boundary_pair_vec resolved_boundaries(const full_hit &,
	                                          const full_hit_list &,
	                                          const crh_segment_spec &);

	seq::seq_seg_opt_vec merge_boundaries(const seq::seq_seg_vec &,
	                                      const seg_boundary_pair_vec &,
	                                      const crh_segment_spec &);

	seq::seq_seg_opt_vec resolve_all_boundaries(const full_hit &,
	                                            const full_hit_list &,
	                                            const crh_segment_spec &);

	std::string get_all_resolved_segments_string(const full_hit &,
	                                             const full_hit_list &,
	                                             const crh_segment_spec &);

	std::string to_json_string_with_compact_fullhits(const full_hit_list &,
	                                                 const crh_segment_spec_opt & = ::std::nullopt,
	                                                 const size_t & = 0);

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FULL_HIT_LIST_FNS_HPP
