/// \file
/// \brief The pdb_base class definitions

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

#include "pdb_base.hpp"

#include "biocore/residue_id.hpp"
#include "structure/geometry/coord.hpp"

using namespace boost::filesystem;
using namespace cath;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

const string pdb_base::PDB_RECORD_STRING_TER ( "TER   " );
constexpr size_t pdb_base::MIN_NUM_PDB_COLS;
constexpr size_t pdb_base::MAX_NUM_PDB_COLS;

/// \brief An NVI pass-through to the virtual do_rotate()
void pdb_base::read_file(const path &arg_filename ///< TODOCUMENT
                         ) {
	do_read_file( arg_filename );
}

/// \brief An NVI pass-through to the virtual do_rotate()
void pdb_base::append_to_file(const path &arg_filename ///< TODOCUMENT
                              ) const {
	do_append_to_file(arg_filename);
}

/// \brief An NVI pass-through to the virtual do_rotate()
void pdb_base::set_chain_label(const chain_label &arg_chain_label ///< TODOCUMENT
                               ) {
	do_set_chain_label( arg_chain_label );
}

/// \brief An NVI pass-through to the virtual do_rotate()
residue_id_vec pdb_base::get_residue_ids_of_first_chain__backbone_unchecked() const {
	return do_get_residue_ids_of_first_chain__backbone_unchecked();
}

/// \brief An NVI pass-through to the virtual do_rotate()
coord pdb_base::get_residue_ca_coord_of_index__backbone_unchecked(const size_t &arg_index ///< TODOCUMENT
                                                                  ) const {
	return do_get_residue_ca_coord_of_index__backbone_unchecked( arg_index );
}

/// \brief An NVI pass-through to the virtual do_get_no_atoms()
size_t pdb_base::get_num_atoms() const {
	return do_get_num_atoms();
}

/// \brief An NVI pass-through to the virtual do_rotate()
void pdb_base::rotate(const rotation &arg_rotation ///< \brief The rotation to apply to this
                      ) {
	do_rotate(arg_rotation);
}

/// \brief An NVI pass-through to the virtual do_add()
void pdb_base::operator+=(const coord &arg_coord ///< \brief The coord to add to this
                          ) {
	do_add(arg_coord);
}

/// \brief An NVI pass-through to the virtual do_subtract()
void pdb_base::operator-=(const coord &arg_coord ///< \brief The coord to subtract from this
                          ) {
	do_subtract(arg_coord);
}

