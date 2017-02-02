/// \file
/// \brief The crh_segment_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_SEGMENT_SPEC_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_SEGMENT_SPEC_H

#include "resolve_hits/trim/trim_spec.hpp"

namespace cath {
	namespace rslv {

		/// \brief Specify how segments should be handled in cath-resolve-hits
		class crh_segment_spec final {
		private:
			/// \brief The specification for trimming hits' segments to allow some overlap
			trim_spec overlap_trim_spec = DEFAULT_OVERLAP_TRIM_SPEC;

			/// \brief The minimum segment length
			residx_t  min_seg_length    = DEFAULT_MIN_SEG_LENGTH;

		public:
			/// \brief The default value for the specification for trimming hits' segments to allow some overlap
			static constexpr trim_spec DEFAULT_OVERLAP_TRIM_SPEC = { 50, 30 };

			/// \brief The default value for the minimum segment length
			static constexpr residx_t  DEFAULT_MIN_SEG_LENGTH    =  7;

			const trim_spec & get_overlap_trim_spec() const;
			const residx_t & get_min_seg_length() const;

			crh_segment_spec & set_overlap_trim_spec(const trim_spec &);
			crh_segment_spec & set_min_seg_length(const residx_t &);
		};

		using crh_segment_spec_opt = boost::optional<crh_segment_spec>;

		crh_segment_spec make_no_action_crh_segment_spec();

		hit_seg_opt apply_spec_to_seg_copy(const hit_seg &,
		                                   const crh_segment_spec_opt &);
	} // namespace rslv
} // namespace cath

#endif
