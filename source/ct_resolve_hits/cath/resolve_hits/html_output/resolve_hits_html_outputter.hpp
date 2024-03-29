/// \file
/// \brief The resolve_hits_html_outputter class header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT_RESOLVE_HITS_HTML_OUTPUTTER_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT_RESOLVE_HITS_HTML_OUTPUTTER_HPP

#include "cath/common/type_aliases.hpp"
#include "cath/resolve_hits/options/spec/crh_filter_spec.hpp"
#include "cath/resolve_hits/options/spec/crh_html_spec.hpp"
#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"

#include <iosfwd>

// clang-format off
namespace cath { class display_colour; }
namespace cath::rslv { class calc_hit_list; }
namespace cath::rslv { class crh_score_spec; }
namespace cath::rslv { class crh_segment_spec; }
namespace cath::rslv { class full_hit; }
namespace cath::rslv { class full_hit_list; }
namespace cath::rslv { class trim_spec; }
// clang-format on

namespace cath::rslv {

	/// \brief The section of the HTML table in which an element appears
	enum class table_section : bool {
		RESULTS, ///< The results section of the table
		INPUTS   ///< The inputs section of the table
	};

	/// \brief The context in which a hit is being rendered into an HTML row
	enum class hit_row_context : char {
		HIGHLIGHT,    ///< An input data row that is being highlighted as containing a hit that appears in the results
		NORMAL,       ///< A normal input data row
		RESULT,       ///< A result row
		RESULT_FULL   ///< A full-result row, combining all the result hits in one row
	};

	std::string row_css_class_and_data_of_hit_row_context(const hit_row_context &);
	std::string first_cell_css_class_of_hit_row_context(const hit_row_context &);

	/// \brief Contain functions for generating HTML to describe cath-resolve hits such as hits, architectures etc
	///
	/// Possible new functions/options
	///  * display up to a maximum of n hits after the last hit that's been chosen in the architecture
	///  * grey out hits that are strictly worse than other hits
	///  * abstract out the style into a CSS
	///  * add JS to remove the regions/length columns and put them into mouseovers
	///  * think about (how/whether) to handle data for multiple query_ids
	class resolve_hits_html_outputter final {
	private:
		static std::string total_score_row(const resscr_t &);
		static std::string markers_row(const size_t &,
		                               const str_opt &,
		                               const table_section &);

		static std::string hit_html(const html_hit &,
		                            const crh_segment_spec &,
		                            const size_t &,
		                            const hit_row_context &);
		static std::string hits_row_html(const html_hit_vec &,
		                                 const crh_segment_spec &,
		                                 const crh_score_spec &,
		                                 const size_t &,
		                                 const hit_row_context &);

	public:
		static std::string html_prefix();
		static std::string html_key();
		static std::string html_suffix();

		static std::string css_string();

		static size_t step_for_length(const size_t &);

		static std::string output_html(const std::string &,
		                               full_hit_list &&,
		                               const crh_score_spec &,
		                               const crh_segment_spec &,
		                               const crh_html_spec & = crh_html_spec{},
		                               const bool & = true,
		                               const crh_filter_spec & = make_accept_all_filter_spec(),
		                               const size_t & = 0);

		static std::string output_html(const std::string &,
		                               const full_hit_list &,
		                               const crh_score_spec &,
		                               const crh_segment_spec &,
		                               const crh_html_spec & = crh_html_spec{},
		                               const bool & = true,
		                               const crh_filter_spec & = make_accept_all_filter_spec(),
		                               const size_t & = 0);

		static std::string output_html(const std::string &,
		                               const calc_hit_list &,
		                               const crh_score_spec &,
		                               const crh_segment_spec &,
		                               const crh_html_spec & = crh_html_spec{},
		                               const bool & = true,
		                               const crh_filter_spec & = make_accept_all_filter_spec(),
		                               const size_t & = 0);

		
	};

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HTML_OUTPUT_RESOLVE_HITS_HTML_OUTPUTTER_HPP
