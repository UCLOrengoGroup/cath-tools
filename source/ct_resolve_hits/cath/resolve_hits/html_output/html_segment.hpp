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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT_HTML_SEGMENT_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT_HTML_SEGMENT_HPP

#include "cath/display_colour/display_colour.hpp"
#include "cath/display_colour/display_colour_type_aliases.hpp"
#include "cath/seq/seq_arrow.hpp"

namespace cath::rslv {

	/// \brief Represent a segment to be rendered in HTML
	struct html_segment final {
	private:

		/// \brief The sides on which a pill should be rounded on
		enum class pill_rounding : char {
			NEITHER,    ///< Round the pill on neither side
			LEFT_ONLY,  ///< Round the pill on the left side
			RIGHT_ONLY, ///< Round the pill on the right side
			BOTH        ///< Round the pill on both sides
		};

		static std::string get_html_string(const seq::seq_arrow &,
		                                   const seq::res_arrow_opt &,
		                                   const std::string &,
		                                   const str_str_pair_vec &,
		                                   const display_colour_opt &,
		                                   const display_colour &,
		                                   const size_t &,
		                                   const pill_rounding & = pill_rounding::BOTH);

		static std::string get_resolve_boundary_html_string(const seq::seq_arrow &,
		                                                    const display_colour &,
		                                                    const size_t &);

	public:
		/// \brief The position of the segment's start
		seq::seq_arrow start;

		/// \brief The position of the segment's trimmed start
		///        (or nullopt if this is so short it isn't counted as a segment)
		seq::res_arrow_opt trimmed_start;

		/// \brief The position of the segment's trimmed stop
		///        (or nullopt if this is so short it isn't counted as a segment)
		seq::res_arrow_opt trimmed_stop;

		/// \brief The position of the segment's stop
		seq::seq_arrow stop;


		/// \brief The position of any resolved boundary start
		///        (ie the meeting point between this segment's start and some other overlapping segment's end)
		seq::res_arrow_opt resolved_start;

		/// \brief The position of any resolved boundary end
		///        (ie the meeting point between this segment's end and some other overlapping segment's start)
		seq::res_arrow_opt resolved_stop;


		/// \brief The colour in which the hit should be rendered
		display_colour colour;

		/// \brief A set of key/value pairs to be inserted as data attributes in the span (keys are prefixed with "data-" if not already)
		str_str_pair_vec data_key_values;

		/// \brief The full length of the sequence on which this hit appears
		size_t full_seq_length;

		[[nodiscard]] std::string get_grey_back_html_string() const;
		[[nodiscard]] std::string get_lightened_back_html_string() const;
		[[nodiscard]] std::string get_full_result_html_string() const;
		[[nodiscard]] std::string get_strong_front_html_string() const;

		[[nodiscard]] str_vec get_all_span_html_strs( const bool & = true ) const;
	};

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT_HTML_SEGMENT_HPP
