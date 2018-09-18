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

#include "crh_score_spec.hpp"

using namespace cath::rslv;

constexpr resscr_t crh_score_spec::DEFAULT_LONG_DOMAINS_PREFERENCE;
constexpr resscr_t crh_score_spec::DEFAULT_HIGH_SCORES_PREFERENCE;
constexpr bool     crh_score_spec::DEFAULT_APPLY_CATH_RULES;
constexpr bool     crh_score_spec::DEFAULT_NAIVE_GREEDY;

/// \brief Ctor
crh_score_spec::crh_score_spec(const bool     &prm_apply_cath_rules,        ///< Whether to apply rules specific to CATH-Gene3D
                               const resscr_t &prm_long_domains_preference, ///< The degree to which long domains are preferred
                               const resscr_t &prm_high_scores_preference,  ///< The degree to which high scores are preferred
                               const bool     &prm_naive_greedy             ///< Whether to use a naive, greedy approach to resolving
                               ) : long_domains_preference { prm_long_domains_preference },
                                   high_scores_preference  { prm_high_scores_preference  },
                                   apply_cath_rules        { prm_apply_cath_rules        },
                                   naive_greedy            { prm_naive_greedy            } {
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

/// \brief Getter for whether to use a naive, greedy approach to resolving
const bool & crh_score_spec::get_naive_greedy() const {
	return naive_greedy;
}

/// \brief Setter for the degree to which long domains are preferred
crh_score_spec & crh_score_spec::set_long_domains_preference(const resscr_t &prm_long_domains_preference ///< The degree to which long domains are preferred
                                                             ) {
	long_domains_preference = prm_long_domains_preference;
	return *this;
}

/// \brief Setter for the degree to which high scores are preferred
crh_score_spec & crh_score_spec::set_high_scores_preference(const resscr_t &prm_high_scores_preference ///< The degree to which high scores are preferred
                                                            ) {
	high_scores_preference = prm_high_scores_preference;
	return *this;
}

/// \brief Setter for whether to apply rules specific to CATH-Gene3D
crh_score_spec & crh_score_spec::set_apply_cath_rules(const bool &prm_apply_cath_rules ///< Whether to apply rules specific to CATH-Gene3D
                                                      ) {
	apply_cath_rules = prm_apply_cath_rules;
	return *this;
}

/// \brief Setter for whether to use a naive, greedy approach to resolving
crh_score_spec & crh_score_spec::set_naive_greedy(const bool &prm_naive_greedy ///< Whether to use a naive, greedy approach to resolving
                                                  ) {
	naive_greedy = prm_naive_greedy;
	return *this;
}

/// \brief Make a neutral crh_score_spec
///
/// \relates crh_score_spec
crh_score_spec cath::rslv::make_neutral_score_spec() {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return crh_score_spec{ false, 0.0, 0.0 };
}
