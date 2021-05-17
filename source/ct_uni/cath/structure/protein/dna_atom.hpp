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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_DNA_ATOM_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_DNA_ATOM_HPP

#include "cath/common/char_arr_type_aliases.hpp"
#include "cath/common/cpp20/arrays_are_equal.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

#include <string>

namespace cath {

	/// \brief Represent the different types of DNA/RNA nucleotide atom that are possible in PDB files
	enum class dna_atom : char {
		A,
		C,
		G,
		I, // eg in 1xnr
		N,
		T, // eg in 3dpv
		U,

		DA,
		DC,
		DG,
		DI, // eg in 3rzl
		DT,
		DU,

		PLUS_A, // eg in 1hp6
		PLUS_C, // eg in 356d
		PLUS_G, // eg in 1gpg
		PLUS_U, // eg in 1hp6

		A_SPACE, // eg in 3zvp
		C_SPACE, // eg in 3zvp
		G_SPACE, // eg in 4a1d
		U_SPACE  // eg in 4a1d
	};

	/// \brief TODOCUMENT
	///
	/// Examples:
	///  * PDBs 1dgo, 1ii1, 1js7, 1qe7, 3ga6 and 4fm0 all contain
	///    ATOM records with each of the residue names: DA, DC, DG, DT and DU
	///  * PDBs 1eg0, 1fjf, 2i82 and 3u5f all contain
	///    ATOM records with each of the residue names: A, C, G, N and U
	constexpr ::std::array DNA_RNA_RESIDUE_NAMES = {
		::std::pair{ make_char_arr( "  A" ), dna_atom::A },
		::std::pair{ make_char_arr( "  C" ), dna_atom::C },
		::std::pair{ make_char_arr( "  G" ), dna_atom::G },
		::std::pair{ make_char_arr( "  I" ), dna_atom::I }, // eg in 1xnr
		::std::pair{ make_char_arr( "  N" ), dna_atom::N },
		::std::pair{ make_char_arr( "  T" ), dna_atom::T }, // eg in 3dpv
		::std::pair{ make_char_arr( "  U" ), dna_atom::U },

		::std::pair{ make_char_arr( " DA" ), dna_atom::DA },
		::std::pair{ make_char_arr( " DC" ), dna_atom::DC },
		::std::pair{ make_char_arr( " DG" ), dna_atom::DG },
		::std::pair{ make_char_arr( " DI" ), dna_atom::DI }, // eg in 3rzl
		::std::pair{ make_char_arr( " DT" ), dna_atom::DT },
		::std::pair{ make_char_arr( " DU" ), dna_atom::DU },

		::std::pair{ make_char_arr( " +A" ), dna_atom::PLUS_A }, // eg in 1hp6
		::std::pair{ make_char_arr( " +C" ), dna_atom::PLUS_C }, // eg in 356d
		::std::pair{ make_char_arr( " +G" ), dna_atom::PLUS_G }, // eg in 1gpg
		::std::pair{ make_char_arr( " +U" ), dna_atom::PLUS_U }, // eg in 1hp6

		::std::pair{ make_char_arr( " A " ), dna_atom::A_SPACE }, // eg in 3zvp
		::std::pair{ make_char_arr( " C " ), dna_atom::C_SPACE }, // eg in 3zvp
		::std::pair{ make_char_arr( " G " ), dna_atom::G_SPACE }, // eg in 4a1d
		::std::pair{ make_char_arr( " U " ), dna_atom::U_SPACE }  // eg in 4a1d
	};

	/// \brief TODOCUMENT
	///
	/// \param prm_letter The 1 letter code
	constexpr ::std::optional<dna_atom> dna_atom_of_code( const char_3_arr &prm_code ) {
		for ( const auto &[ code, dna ] : DNA_RNA_RESIDUE_NAMES ) {
			if ( common::arrays_are_equal( code, prm_code ) ) {
				return { dna };
			}
		}
		return ::std::nullopt;
	}

	/// \brief TODOCUMENT
	///
	/// \param prm_dna_atom The dna_atom value
	constexpr char_3_arr to_three_char_arr( const dna_atom &prm_dna_atom ) {
		for ( const auto &[ code, dna ] : DNA_RNA_RESIDUE_NAMES ) {
			if ( dna == prm_dna_atom ) {
				return code;
			}
		}
		/// TODO: Come GCC >= 10, remove this silly dance to appeas it about the unconditionally non-constexpr throwing after here
		if ( DNA_RNA_RESIDUE_NAMES.front().second != prm_dna_atom ) {
			return DNA_RNA_RESIDUE_NAMES.front().first;
		}
		else {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unhandled dna_atom"));
		}
	}

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_DNA_ATOM_HPP
