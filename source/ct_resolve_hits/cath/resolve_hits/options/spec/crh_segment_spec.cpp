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

using namespace ::cath::rslv;
using namespace ::cath::seq;

using ::std::make_optional;
using ::std::nullopt;

constexpr trim_spec crh_segment_spec::DEFAULT_OVERLAP_TRIM_SPEC;
constexpr residx_t  crh_segment_spec::DEFAULT_MIN_SEG_LENGTH;

/// \brief Getter for the specification for trimming hits' segments to allow some overlap
const trim_spec & crh_segment_spec::get_overlap_trim_spec() const {
	return overlap_trim_spec;
}

/// \brief Getter for the minimum segment length
const residx_t & crh_segment_spec::get_min_seg_length() const {
	return min_seg_length;
}

/// \brief Setter for the specification for trimming hits' segments to allow some overlap
crh_segment_spec & crh_segment_spec::set_overlap_trim_spec(const trim_spec &prm_overlap_trim_spec ///< The specification for trimming hits' segments to allow some overlap
                                                           ) {
	overlap_trim_spec = prm_overlap_trim_spec;
	return *this;
}

/// \brief Setter for the minimum segment length
crh_segment_spec & crh_segment_spec::set_min_seg_length(const residx_t &prm_min_seg_length ///< The minimum segment length
                                                        ) {
	min_seg_length = prm_min_seg_length;
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

/// \brief Apply the specified crh_segment_spec to a copy of the specified seq_seg, returning
///        a trimmed version if the segment's length meets the min-seg-length, or nullopt otherwise
///
/// \relates crh_segment_spec
seq_seg_opt cath::rslv::apply_spec_to_seg_copy(const seq_seg              &prm_seg,         ///< The seq_seg to copy
                                               const crh_segment_spec_opt &prm_segment_spec ///< The crh_segment_spec to apply
                                               ) {
	return 
		( prm_segment_spec && get_length( prm_seg ) >= prm_segment_spec->get_min_seg_length() )
			? make_optional( trim_seq_seg_copy( prm_seg, prm_segment_spec->get_overlap_trim_spec() ) )
			: nullopt;
}

/// \brief Get an optional trim_spec of the specified optional crh_segment_spec
///        (where the result is nullopt iff the input is nullopt)
///
/// \relates crh_segment_spec
trim_spec_opt cath::rslv::get_trim_spec_opt(const crh_segment_spec_opt &prm_crh_segment_spec ///< The optional crh_segment_spec from which to make an optional trim_spec
                                            ) {
	return prm_crh_segment_spec
		? make_optional( prm_crh_segment_spec->get_overlap_trim_spec() )
		: nullopt;
}
