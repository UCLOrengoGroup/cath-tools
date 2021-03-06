/// \file
/// \brief The wolf_file class definitions

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

#include "wolf_file.hpp"

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::std;

/// \brief Ctor for wolf_file
wolf_file::wolf_file(residue_vec prm_wolf_residues
                     ) : wolf_residues{ std::move( prm_wolf_residues ) } {
}

/// \brief Return the number of residues
wolf_file::size_type wolf_file::get_num_residues() const {
	return wolf_residues.size();
}

/// \brief Return the residue at the specified index
const residue & wolf_file::get_residue_of_index(const size_type &prm_index ///< The index of the residue to return
                                                ) const {
	return wolf_residues[prm_index];
}

/// \brief Convert a wolf_file to a protein
///
/// \relates wolf_file
///
/// \relates protein
protein cath::file::protein_from_wolf(const wolf_file &prm_wolf_file ///< The wolf_file object for a given structure
	                                  ) {
	// Copy the residues from the wolf_file
	const wolf_file::size_type num_wolf_residues = prm_wolf_file.get_num_residues();
	residue_vec wolf_residues;
	wolf_residues.reserve(num_wolf_residues);
	for (const wolf_file::size_type &wolf_residue_ctr : indices( num_wolf_residues ) ) {
		wolf_residues.push_back(prm_wolf_file.get_residue_of_index(wolf_residue_ctr));
	}

	// Put the residues in a new protein object and return it
	protein new_protein;
	new_protein.set_residues(wolf_residues);
	return new_protein;
}

