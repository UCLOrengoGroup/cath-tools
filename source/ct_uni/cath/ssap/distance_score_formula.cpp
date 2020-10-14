/// \file
/// \brief The distance_score_formula definitions

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

#include "distance_score_formula.hpp"

#include <iostream>

using namespace ::cath;
using namespace ::std;

/// \brief TODOCUMENT
map<distance_score_formula, string> name_of_distance_score_formula::get() {
	return {
		{ distance_score_formula::FROM_SSAP_PAPER,       "distance_score_formula_from_ssap_paper"       },
		{ distance_score_formula::USED_IN_PREVIOUS_CODE, "distance_score_formula_used_in_previous_code" },
		{ distance_score_formula::SIMPLIFIED,            "distance_score_formula_simplified"            }
	};
}

/// \brief TODOCUMENT
///
/// \relates distance_score_formula
ostream & cath::operator<<(ostream                      &prm_os,                    ///< TODOCUMENT
                           const distance_score_formula &prm_distance_score_formula ///< TODOCUMENT
						   ) {
	prm_os << name_of_distance_score_formula::get().at( prm_distance_score_formula );
	return prm_os;
}
