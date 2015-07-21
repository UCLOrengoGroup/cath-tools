/// \file
/// \brief The selection_policy_acquirer class definitions

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

#include "selection_policy_acquirer.h"

#include "alignment/alignment_coord_extractor.h"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_list.h"
#include "structure/geometry/coord.h"
#include "structure/geometry/coord_list.h"
#include "structure/protein/protein.h"
#include "structure/protein/residue.h"

using namespace cath::align;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::opts;
using namespace std;

/// \brief TODOCUMENT
const common_residue_selection_policy & selection_policy_acquirer::get_comm_res_seln_pol_ptr_cref() const {
	return *comm_res_seln_pol_ptr;
}

/// \brief TODOCUMENT
const common_atom_selection_policy & selection_policy_acquirer::get_comm_atom_seln_pol_ptr_cref() const {
	return *comm_atom_seln_pol_ptr;
}

/// \brief Ctor for selection_policy_acquirer
selection_policy_acquirer::selection_policy_acquirer(const common_residue_selection_policy &arg_comm_res_seln_pol, ///< TODOCUMENT
                                                     const common_atom_selection_policy    &arg_comm_atom_seln_pol ///< TODOCUMENT
                                                     ) : comm_res_seln_pol_ptr(   arg_comm_res_seln_pol.clone() ),
                                                         comm_atom_seln_pol_ptr( arg_comm_atom_seln_pol.clone() ) {
}

/// \brief TODOCUMENT
pair<coord_list, coord_list> selection_policy_acquirer::get_common_coords(const alignment &arg_alignment,   ///< TODOCUMENT
                                                                          const pdb       &arg_pdb_1,       ///< TODOCUMENT
                                                                          const pdb       &arg_pdb_2,       ///< TODOCUMENT
//                                                                          const str_vec              &/*arg_names*/, ///< TODOCUMENT
                                                                          const size_t    &arg_index_1,     ///< TODOCUMENT
                                                                          const size_t    &arg_index_2      ///< TODOCUMENT
                                                                          ) const {
	return alignment_coord_extractor::get_common_coords(
		arg_alignment,
		arg_pdb_1,
		arg_pdb_2,
		get_comm_res_seln_pol_ptr_cref(),
		get_comm_atom_seln_pol_ptr_cref(),
		arg_index_1,
		arg_index_2
	);
}

/// \brief TODOCUMENT
string selection_policy_acquirer::get_descriptive_name() const {
	return get_comm_res_seln_pol_ptr_cref().get_descriptive_name()
	       + "."
	       + get_comm_atom_seln_pol_ptr_cref().get_descriptive_name();
}

/// \brief TODOCUMENT
coord_list_coord_list_pair cath::opts::get_common_coords(const selection_policy_acquirer &arg_seln_pol_acq, ///< TODOCUMENT
                                                         const alignment                 &arg_alignment,    ///< TODOCUMENT
                                                         const pdb_list                  &arg_pdbs,         ///< TODOCUMENT
                                                         // const str_vec S                  &arg_names,        ///< TODOCUMENT
                                                         const size_t                    &arg_index_1,      ///< TODOCUMENT
                                                         const size_t                    &arg_index_2       ///< TODOCUMENT
                                                         ) {
	return arg_seln_pol_acq.get_common_coords(
		arg_alignment,
		arg_pdbs[ arg_index_1 ],
		arg_pdbs[ arg_index_2 ],
		arg_index_1,
		arg_index_2
	);
}

