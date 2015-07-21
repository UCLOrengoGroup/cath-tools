/// \file
/// \brief The tally_residue_names class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef TALLY_RESIDUE_NAMES_H_INCLUDED
#define TALLY_RESIDUE_NAMES_H_INCLUDED

#include "common/type_aliases.h"
#include "structure/structure_type_aliases.h"

namespace cath {
	namespace file {

		size_size_pair_vec tally_residue_names(const residue_name_vec &,
		                                       const residue_name_vec &,
		                                       const bool &,
		                                       const bool &arg_permit_tail_break_without_null_residue = true);

		size_size_pair_vec tally_residue_names_str(const residue_name_vec &,
		                                           const str_vec &,
		                                           const bool &,
		                                           const bool &arg_permit_tail_break_without_null_residue = true);
	}
}

#endif
