/// \file
/// \brief The quad_find_action_check class definitions

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

#include "quad_find_action_check.hpp"

using namespace cath::index;
using namespace cath::index::detail;
// using namespace std;

/// \brief Ctor for quad_find_action_check
quad_find_action_check::quad_find_action_check(const protein             &arg_protein_a, ///< The first  protein being compared
                                               const protein             &arg_protein_b, ///< The second protein being compared
                                               const vcie_match_criteria &arg_criteria   ///< The criteria for which entries should be compared
                                               ) : protein_a    ( arg_protein_a ),
                                                   protein_b    ( arg_protein_b ),
                                                   the_criteria ( arg_criteria  ) {
}

/// \brief Getter for the total score
const double & quad_find_action_check::get_total_score() const {
	return total_score;
}
