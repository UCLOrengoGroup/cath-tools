/// \file
/// \brief The tally_residue_ids class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_DSSP_WOLF_TALLY_RESIDUE_NAMES_H
#define _CATH_TOOLS_SOURCE_FILE_DSSP_WOLF_TALLY_RESIDUE_NAMES_H

#include "common/type_aliases.hpp"
#include "structure/structure_type_aliases.hpp"

namespace cath {
	namespace file {

		size_size_pair_vec tally_residue_ids(const residue_id_vec &,
		                                     const residue_id_vec &,
		                                     const bool &,
		                                     const bool & = true,
		                                     const size_set & = size_set{} );

	} // namespace file
} // namespace cath

#endif
