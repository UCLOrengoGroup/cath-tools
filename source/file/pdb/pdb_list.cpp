/// \file
/// \brief The pdb_list class definitions

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

#include "pdb_list.h"

#include "common/c++14/cbegin_cend.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_list.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

/// \brief TODOCUMENT
void pdb_list::push_back(const pdb &arg_pdb ///< TODOCUMENT
                         ) {
	pdbs.push_back(arg_pdb);
}

/// \brief TODOCUMENT
void pdb_list::reserve(const size_t &arg_size ///< TODOCUMENT
                       ) {
	pdbs.reserve(arg_size);
}

/// \brief TODOCUMENT
size_t pdb_list::size() const {
	return pdbs.size();
}

/// \brief TODOCUMENT
bool pdb_list::empty() const {
	return pdbs.empty();
}

/// \brief TODOCUMENT
pdb & pdb_list::operator[](const size_t &arg_index ///< TODOCUMENT
                           ) {
	return pdbs[arg_index];
}

/// \brief TODOCUMENT
const pdb & pdb_list::operator[](const size_t &arg_index ///< TODOCUMENT
                                 ) const {
	return pdbs[arg_index];
}

///// \brief TODOCUMENT
//pdb_list::iterator pdb_list::begin() {
//	return std::begin( pdbs );
//}
///// \brief TODOCUMENT
//pdb_list::iterator pdb_list::end() {
//	return std::end( pdbs );
//}
/// \brief TODOCUMENT
pdb_list::const_iterator pdb_list::begin() const {
	return common::cbegin( pdbs );
}
/// \brief TODOCUMENT
pdb_list::const_iterator pdb_list::end() const {
	return common::cend( pdbs );
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
pdb_list cath::file::make_pdb_list(const pdb_vec &arg_pdbs ///< TODOCUMENT
                                   ) {
	pdb_list new_pdb_list;
	new_pdb_list.reserve( arg_pdbs.size() );
	for (const pdb &the_pdb : arg_pdbs) {
		new_pdb_list.push_back( the_pdb );
	}
	return new_pdb_list;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
pdb_list cath::file::pdb_list_of_backbone_complete_subset_pdbs(const pdb_list &arg_pdb_list, ///< TODOCUMENT
                                                               ostream        &arg_ostream   ///< TODOCUMENT
                                                               ) {
	pdb_list new_pdb_list;
	new_pdb_list.reserve( arg_pdb_list.size() );
	for (const pdb &the_pdb : arg_pdb_list) {
		new_pdb_list.push_back( backbone_complete_subset_of_pdb( the_pdb, arg_ostream ) );
	}
	return new_pdb_list;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
///
/// \relates protein_list
protein_list cath::file::build_protein_list_of_pdb_list(const pdb_list &arg_pdb_list ///< TODOCUMENT
                                                        ) {
	protein_list new_protein_list;
	new_protein_list.reserve(arg_pdb_list.size());
	for (const pdb &the_pdb : arg_pdb_list) {
		new_protein_list.push_back( build_protein_of_pdb( the_pdb ) );
	}
	return new_protein_list;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
///
/// \relates protein_list
protein_list cath::file::build_protein_list_of_pdb_list_and_names(const pdb_list &arg_pdb_list, ///< TODOCUMENT
                                                                  const str_vec  &arg_names     ///< TODOCUMENT
                                                                  ) {
	const size_t num_names = arg_names.size();
	if ( arg_pdb_list.size() != num_names ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to make proteins from pdb_list and names because the numbers mismatch"));
	}
	protein_list new_proteins = build_protein_list_of_pdb_list( arg_pdb_list );
	for (size_t protein_ctr = 0; protein_ctr < num_names; ++protein_ctr) {
		new_proteins[ protein_ctr ].set_title( arg_names[ protein_ctr ] );
	}
	return new_proteins;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
amino_acid_vec_vec cath::file::get_amino_acid_lists(const pdb_list &arg_pdb_list ///< TODOCUMENT
                                                    ) {
	amino_acid_vec_vec amino_acid_lists;
	amino_acid_lists.reserve( arg_pdb_list.size() );
	for (const pdb &the_pdb : arg_pdb_list) {
		amino_acid_lists.push_back( get_amino_acid_list( the_pdb ) );
	}
	return amino_acid_lists;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
residue_name_vec_vec cath::file::get_residue_names_of_first_chains__backbone_unchecked(const pdb_list &arg_pdb_list ///< TODOCUMENT
                                                                                       ) {
	const size_t num_pdbs = arg_pdb_list.size();
	residue_name_vec_vec residue_names;
	residue_names.reserve( num_pdbs );
	for (size_t pdb_ctr = 0; pdb_ctr < num_pdbs; ++pdb_ctr) {
		residue_names.push_back( arg_pdb_list[ pdb_ctr ].get_residue_names_of_first_chain__backbone_unchecked() );
	}
	return residue_names;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
residue_name_vec_vec cath::file::get_backbone_complete_residue_names_of_first_chains(const pdb_list &arg_pdb_list ///< TODOCUMENT
                                                                                     ) {
	const size_t num_pdbs = arg_pdb_list.size();
	residue_name_vec_vec residue_names;
	residue_names.reserve( num_pdbs );
	for (size_t pdb_ctr = 0; pdb_ctr < num_pdbs; ++pdb_ctr) {
		residue_names.push_back( arg_pdb_list[ pdb_ctr ].get_backbone_complete_residue_names_of_first_chain() );
	}
	return residue_names;
}