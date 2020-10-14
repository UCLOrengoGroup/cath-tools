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

#include "cath/chopping/region/region.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_list.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/superposition/superposition.hpp"
#include "cath/superposition/superposition_context.hpp"

#include <string>

using namespace cath::opts;
using namespace cath::sup;

using std::ostream;

constexpr double superposition_acquirer::TOLERANCE_FOR_EQUAL_RMSDS;

/// \brief TODOCUMENT
superposition_context superposition_acquirer::get_superposition(ostream &prm_stderr ///< TODOCUMENT
                                                                ) const {
	return do_get_superposition( prm_stderr );
}
