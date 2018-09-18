/// \file
/// \brief The score_classn_value_better_value class definitions

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

#include "score_classn_value_better_value.hpp"

#include "score/score_classification/score_classn_value.hpp"

#include <limits>

using namespace cath::score;
using namespace std;

/// \brief Ctor for that allows the user to specify higher_is_better
score_classn_value_better_value::score_classn_value_better_value(const bool &prm_higher_is_better ///< Whether a higher score_value is "better" and should be treated as less-than
                                                                 ) : higher_is_better( prm_higher_is_better ) {
}

/// \brief Getter for higher_is_better (whether a higher score_value is "better" and should be treated as less-than)
const bool & score_classn_value_better_value::get_higher_is_better() const {
	return higher_is_better;
}

/// \brief Function operator to return whether the first score_classn_value has a "better" score_value than the second
bool score_classn_value_better_value::operator()(const score_classn_value &prm_score_classn_value_a, ///< The first  score_classn_value to compare
                                                 const score_classn_value &prm_score_classn_value_b  ///< The second score_classn_value to compare
                                                 ) {
	// Grab const-references to the two score values
	const double &score_value_a = prm_score_classn_value_a.get_score_value();
	const double &score_value_b = prm_score_classn_value_b.get_score_value();

	// Return whether the first score_value_a is better than the first
	return get_higher_is_better() ? ( score_value_a > score_value_b)
	                              : ( score_value_a < score_value_b);
}

/// \brief TODOCUMENT
double cath::score::get_worst_possible_score(const score_classn_value_better_value &prm_score_classn_value_better_value ///< TODOCUMENT
                                             ) {
	return get_worst_possible_score( prm_score_classn_value_better_value.get_higher_is_better() );
}

/// \brief TODOCUMENT
double cath::score::get_worst_possible_score(const bool &prm_higher_is_better ///< TODOCUMENT
                                             ) {
	return prm_higher_is_better ? numeric_limits<double>::lowest()
	                            : numeric_limits<double>::max();
}
