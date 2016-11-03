/// \file
/// \brief The html_segment class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_HTML_OUTPUT_HTML_SEGMENT_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_HTML_OUTPUT_HTML_SEGMENT_H

#include <boost/optional.hpp>

#include "display_colour/display_colour.h"
#include "display_colour/display_colour_type_aliases.h"
#include "resolve_hits/res_arrow.h"

namespace cath {
	namespace rslv {

		/// \brief Represent a segment to be rendered in HTML
		struct html_segment final {
		private:

			/// \brief The sides on which a pill should be rounded on
			enum class pill_rounding {
				NEITHER,    ///< Round the pill on neither side
				LEFT_ONLY,  ///< Round the pill on the left side
				RIGHT_ONLY, ///< Round the pill on the right side
				BOTH        ///< Round the pill on both sides
			};

			static std::string get_html_string(const res_arrow &,
			                                   const res_arrow_opt &,
			                                   const str_vec &,
			                                   const display_colour_opt &,
			                                   const display_colour &,
			                                   const size_t &,
			                                   const pill_rounding & = pill_rounding::BOTH);

			static std::string get_resolve_boundary_html_string(const res_arrow &,
			                                                    const display_colour &,
			                                                    const size_t &);

		public:
			/// \brief The position of the segment's start
			res_arrow start;

			/// \brief The position of the segment's trimmed start
			///        (or none if this is so short it isn't counted as a segment)
			res_arrow_opt trimmed_start;

			/// \brief The position of the segment's trimmed stop
			///        (or none if this is so short it isn't counted as a segment)
			res_arrow_opt trimmed_stop;

			/// \brief The position of the segment's stop
			res_arrow stop;


			/// \brief The position of any resolved boundary start
			///        (ie the meeting point between this segment's start and some other overlapping segment's end)
			res_arrow_opt resolved_start;

			/// \brief The position of any resolved boundary end
			///        (ie the meeting point between this segment's end and some other overlapping segment's start)
			res_arrow_opt resolved_stop;


			/// \brief The colour in which the hit should be rendered
			display_colour colour;

			/// \brief The index of the hit
			///        (could potentially used for Javascript to identify hits)
			size_t hit_idx;

			/// \brief The full length of the sequence on which this hit appears
			size_t full_seq_length;

			std::string get_grey_back_html_string() const;
			std::string get_lightened_back_html_string() const;
			std::string get_strong_front_html_string() const;

			str_vec get_all_span_html_strs() const;
		};

	} // namespace rslv
} // namespace cath

#endif
