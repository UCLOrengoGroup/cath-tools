/// \file
/// \brief The dna_atom enum header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_DNA_ATOM_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_DNA_ATOM_H

namespace cath {

	/// \brief Represent the different types of DNA/RNA nucleotide atom that are possible in PDB files
	enum class dna_atom : char {
		DA, DC, DG, DT, DU,
		 A,  C,  G,  N,  U
	};

	std::string to_three_char_str(const dna_atom &);

} // namespace cath

#endif
