/// \file
/// \brief The crh_segment_spec class definitions

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

#include "crh_segment_spec.hpp"

using namespace cath::rslv;

constexpr trim_spec           crh_segment_spec::DEFAULT_OVERLAP_TRIM_SPEC;
constexpr residx_t            crh_segment_spec::DEFAULT_MIN_SEG_LENGTH;

/// \brief Getter for the specification for trimming hits' segments to allow some overlap
const trim_spec & crh_segment_spec::get_overlap_trim_spec() const {
	return overlap_trim_spec;
}

/// \brief Getter for the minimum segment length
const residx_t & crh_segment_spec::get_min_seg_length() const {
	return min_seg_length;
}

/// \brief Setter for the specification for trimming hits' segments to allow some overlap
crh_segment_spec & crh_segment_spec::set_overlap_trim_spec(const trim_spec &arg_overlap_trim_spec ///< The specification for trimming hits' segments to allow some overlap
                                                           ) {
	overlap_trim_spec = arg_overlap_trim_spec;
	return *this;
}

/// \brief Setter for the minimum segment length
crh_segment_spec & crh_segment_spec::set_min_seg_length(const residx_t &arg_min_seg_length ///< The minimum segment length
                                                        ) {
	min_seg_length = arg_min_seg_length;
	return *this;
}

/// \brief Make a crh_segment_spec that does nothing to segments
///
/// \relates crh_segment_spec
crh_segment_spec cath::rslv::make_no_action_crh_segment_spec() {
	return crh_segment_spec{}
		.set_overlap_trim_spec( make_no_trim_trim_spec() )
		.set_min_seg_length   ( 1                        );
}

/// \brief Apply the specified crh_segment_spec to a copy of the specified hit_seg, returning
///        a trimmed version if the segment's length meets the min-seg-length, or none otherwise
///
/// \relates crh_segment_spec
hit_seg_opt cath::rslv::apply_spec_to_seg_copy(const hit_seg          &arg_seg,         ///< The hit_seg to copy
                                               const crh_segment_spec &arg_segment_spec ///< The crh_segment_spec to apply
                                               ) {
	return make_optional(
		get_length( arg_seg ) >= arg_segment_spec.get_min_seg_length(),
		trim_hit_seg_copy( arg_seg, arg_segment_spec.get_overlap_trim_spec() )
	);
}