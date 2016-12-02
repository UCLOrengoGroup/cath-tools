/// \file
/// \brief The superposition_acquirer class definitions

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

#include "superposition_acquirer.hpp"

//#include <boost/test/floating_point_comparison.hpp>

//#include "alignment/alignment_coord_extractor.hpp"
//#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.hpp"
//#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
//#include "exception/invalid_argument_exception.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
//#include "structure/geometry/coord_list.hpp"
#include "superposition/superposition.hpp"
#include "superposition/superposition_context.hpp"

//#include <iostream> // **** TEMPORARY ****
#include <string>

//using namespace boost;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

constexpr double superposition_acquirer::PERCENT_TOLERANCE_FOR_EQUAL_RMSDS;

/// \brief TODOCUMENT
superposition_context superposition_acquirer::get_superposition(ostream &arg_stderr ///< TODOCUMENT
                                                                ) const {
	return do_get_superposition( arg_stderr );
}
