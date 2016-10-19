/// \file
/// \brief The crh_score_spec class definitions

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

#include "crh_score_spec.h"

using namespace cath::rslv;

constexpr resscr_t crh_score_spec::DEFAULT_LONG_DOMAINS_PREFERENCE;
constexpr resscr_t crh_score_spec::DEFAULT_HIGH_SCORES_PREFERENCE;
constexpr bool     crh_score_spec::DEFAULT_APPLY_CATH_RULES;

/// \brief Ctor
crh_score_spec::crh_score_spec(const bool     &arg_apply_cath_rules,        ///< Whether to apply rules specific to CATH-Gene3D
                               const resscr_t &arg_long_domains_preference, ///< The degree to which long domains are preferred
                               const resscr_t &arg_high_scores_preference   ///< The degree to which high scores are preferred
                               ) : long_domains_preference { arg_long_domains_preference },
                                   high_scores_preference  { arg_high_scores_preference  },
                                   apply_cath_rules        { arg_apply_cath_rules        } {
}

/// \brief Getter for the degree to which long domains are preferred
const resscr_t & crh_score_spec::get_long_domains_preference() const {
	return long_domains_preference;
}

/// \brief Getter for the degree to which high scores are preferred
const resscr_t & crh_score_spec::get_high_scores_preference() const {
	return high_scores_preference;
}


/// \brief Getter for whether to apply rules specific to CATH-Gene3D
const bool & crh_score_spec::get_apply_cath_rules() const {
	return apply_cath_rules;
}

/// \brief Setter for the degree to which long domains are preferred
crh_score_spec & crh_score_spec::set_long_domains_preference(const resscr_t &arg_long_domains_preference ///< The degree to which long domains are preferred
                                                             ) {
	long_domains_preference = arg_long_domains_preference;
	return *this;
}

/// \brief Setter for the degree to which high scores are preferred
crh_score_spec & crh_score_spec::set_high_scores_preference(const resscr_t &arg_high_scores_preference ///< The degree to which high scores are preferred
                                                            ) {
	high_scores_preference = arg_high_scores_preference;
	return *this;
}

/// \brief Setter for whether to apply rules specific to CATH-Gene3D
crh_score_spec & crh_score_spec::set_apply_cath_rules(const bool &arg_apply_cath_rules ///< Whether to apply rules specific to CATH-Gene3D
                                                      ) {
	apply_cath_rules = arg_apply_cath_rules;
	return *this;
}

/// \brief Make a neutral crh_score_spec
///
/// \relates crh_score_spec
crh_score_spec cath::rslv::make_neutral_score_spec() {
	return { false, 0.0, 0.0 };
}
