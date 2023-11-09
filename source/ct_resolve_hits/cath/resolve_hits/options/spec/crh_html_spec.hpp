/// \file
/// \brief The crh_html_spec class header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_CRH_HTML_SPEC_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_CRH_HTML_SPEC_HPP

#include "cath/common/path_type_aliases.hpp"
#include "cath/resolve_hits/file/hits_input_format_tag.hpp"
#include "cath/resolve_hits/options/spec/crh_html_spec.hpp"
#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"

namespace cath::rslv {

	/// \brief Specify the HTML rendering for cath-resolve-hits
	class crh_html_spec final {
	private:
		/// \brief Whether to restrict HTML output to the contents of the body tag
		bool   restrict_html_within_body = DEFAULT_RESTRICT_HTML_WITHIN_BODY;

		/// \brief The maximum number of non-solution hits to display in the HTML
		size_t max_num_non_soln_hits     = DEFAULT_MAX_NUM_NON_SOLN_HITS;

		/// \brief Whether to exclude hits rejected by the score filters from the HTML
		bool   exclude_rejected_hits     = DEFAULT_EXCLUDE_REJECTED_HITS;

	public:
		/// \brief The default value for whether to restrict HTML output to the contents of the body tag
		static constexpr bool   DEFAULT_RESTRICT_HTML_WITHIN_BODY = false;

		/// \brief The default value for the maximum number of non-solution hits to display in the HTML
		static constexpr size_t DEFAULT_MAX_NUM_NON_SOLN_HITS     = 80;

		/// \brief The default value for whether to exclude hits rejected by the score filters from the HTML
		static constexpr bool   DEFAULT_EXCLUDE_REJECTED_HITS     = false;

		[[nodiscard]] const bool &  get_restrict_html_within_body() const;
		[[nodiscard]] const size_t &get_max_num_non_soln_hits() const;
		[[nodiscard]] const bool &  get_exclude_rejected_hits() const;

		crh_html_spec & set_restrict_html_within_body(const bool &);
		crh_html_spec & set_max_num_non_soln_hits(const size_t &);
		crh_html_spec & set_exclude_rejected_hits(const bool &);
	};

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_CRH_HTML_SPEC_HPP
