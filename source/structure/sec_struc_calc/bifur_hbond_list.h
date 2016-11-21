/// \file
/// \brief The bifur_hbond_list class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_SEC_STRUC_CALC_BIFUR_HBOND_LIST_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_SEC_STRUC_CALC_BIFUR_HBOND_LIST_H

// #include "file/pdb/pdb.h"
// #include "file/pdb/pdb_residue.h"
// #include "structure/geometry/coord.h"
// #include "structure/protein/amino_acid.h"

#include <utility>

namespace cath {
	namespace sec {

		using hbond_partner_t = unsigned int;

		/// \brief TODOCUMENT
		struct dssp_hbond_half final {
			hbond_partner_t index;
			double energy;
		};

		/// \brief TODOCUMENT
		using dssp_hbond_half_pair = std::pair<dssp_hbond_half, dssp_hbond_half>;

		// struct bifur_dssp_hbond_list final {
		// 	dssp_hbond_half_pair
		// };

		// struct bifur_dssp_hbond_list final {
		// 	dssp_hbond_half_pair
		// };
		

	}
} // namespace cath

#endif
