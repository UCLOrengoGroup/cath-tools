/// \file
/// \brief The crh_html_spec class definitions

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

#include "crh_html_spec.hpp"

using namespace cath::rslv;

constexpr bool   crh_html_spec::DEFAULT_RESTRICT_HTML_WITHIN_BODY;
constexpr size_t crh_html_spec::DEFAULT_MAX_NUM_NON_SOLN_HITS;
constexpr bool   crh_html_spec::DEFAULT_EXCLUDE_REJECTED_HITS;

/// \brief Getter for whether to restrict HTML output to the contents of the body tag
const bool & crh_html_spec::get_restrict_html_within_body() const {
	return restrict_html_within_body;
}

/// \brief Getter for the maximum number of non-solution hits to display in the HTML
const size_t & crh_html_spec::get_max_num_non_soln_hits() const {
	return max_num_non_soln_hits;
}

/// \brief Getter for whether to exclude hits rejected by the score filters from the HTML
const bool & crh_html_spec::get_exclude_rejected_hits() const {
	return exclude_rejected_hits;
}

/// \brief Setter for whether to restrict HTML output to the contents of the body tag
crh_html_spec & crh_html_spec::set_restrict_html_within_body(const bool &prm_restrict_html_within_body ///< Whether to restrict HTML output to the contents of the body tag
                                                             ) {
	restrict_html_within_body = prm_restrict_html_within_body;
	return *this;
}

/// \brief Setter for the maximum number of non-solution hits to display in the HTML
crh_html_spec & crh_html_spec::set_max_num_non_soln_hits(const size_t &prm_max_num_non_soln_hits ///< The maximum number of non-solution hits to display in the HTML
                                                         ) {
	max_num_non_soln_hits = prm_max_num_non_soln_hits;
	return *this;
}

/// \brief Setter for whether to exclude hits rejected by the score filters from the HTML
crh_html_spec & crh_html_spec::set_exclude_rejected_hits(const bool &prm_exclude_rejected_hits ///< Whether to exclude hits rejected by the score filters from the HTML
                                                         ) {
	exclude_rejected_hits = prm_exclude_rejected_hits;
	return *this;
}
