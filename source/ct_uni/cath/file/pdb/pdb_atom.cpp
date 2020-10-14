/// \file
/// \brief The pdb_atom class definitions

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

#include "pdb_atom.hpp"

#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "cath/structure/geometry/rotation.hpp"

#include <iomanip>
#include <sstream>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::geom;
using namespace ::std;

using ::boost::format;
using ::boost::io::ios_flags_saver;

constexpr size_t pdb_atom::MIN_NUM_PDB_COLS;
constexpr size_t pdb_atom::MAX_NUM_PDB_COLS;

/// \brief TODOCUMENT
void pdb_atom::rotate(const rotation &prm_rotation ///< TODOCUMENT
                      ) {
	cath::geom::rotate( prm_rotation, atom_coord );
}

/// \brief TODOCUMENT
void pdb_atom::operator+=(const coord &prm_coord ///< TODOCUMENT
                          ) {
	atom_coord += prm_coord;
}

/// \brief TODOCUMENT
void pdb_atom::operator-=(const coord &prm_coord ///< TODOCUMENT
                          ) {
	atom_coord -= prm_coord;
}

/// \brief Get the one letter amino acid code (eg 'S') from a pdb_atom
char cath::file::get_amino_acid_letter_tolerantly(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
                                                  ) {
	return prm_pdb_atom.get_amino_acid().get_letter_tolerantly();
}

/// \brief Get the three letter amino acid code (eg "SER") char_3_arr from a pdb_atom
char_3_arr cath::file::get_amino_acid_code(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
                                           ) {
	return prm_pdb_atom.get_amino_acid().get_code();
}

/// \brief Get the three letter amino acid code (eg "SER") string from a pdb_atom
string cath::file::get_amino_acid_code_string(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
                                       ) {
	return get_code_string( prm_pdb_atom.get_amino_acid() );
}

/// \brief Get the three letter amino acid code (eg "SER") from a pdb_atom
string cath::file::get_amino_acid_name(const pdb_atom &prm_pdb_atom ///< The pdb_atom to query
                                       ) {
	return prm_pdb_atom.get_amino_acid().get_name();
}

/// \brief TODOCUMENT
///
/// \relates pdb_atom
ostream & cath::file::write_pdb_file_entry(ostream          &prm_os,      ///< TODOCUMENT
                                           const residue_id &prm_res_id,  ///< TODOCUMENT
                                           const pdb_atom   &prm_pdb_atom ///< TODOCUMENT
                                           ) {
	// Sanity check the inputs
	if ( prm_res_id.get_residue_name().is_null() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Empty residue_name in cath::write_pdb_file_entry()"));
	}

	// Save the state of the ostream's flags (and reset them in this guard's dtor)
	const ios_flags_saver stream_flags_guard{ prm_os };

	// Grab the necessary data
	const coord  atom_coord                        = prm_pdb_atom.get_coord();
	const string residue_name_with_insert_or_space = make_residue_name_string_with_insert_or_space(
		prm_res_id.get_residue_name()
	);

	                                                                                       // Comments with PDB format documentation
	                                                                                       // (http://www.wwpdb.org/documentation/format33/sect9.html#ATOM)
	prm_os << left;
	prm_os << setw( 6 ) << prm_pdb_atom.get_record_type();                                 //  1 -  6        Record name   "ATOM  " or "HETATM"
	prm_os << right;
	prm_os << setw( 5 ) << prm_pdb_atom.get_atom_serial();                                 //  7 - 11        Integer       serial       Atom  serial number.
	prm_os << " ";
	prm_os << setw( 4 ) << get_element_type_untrimmed_str_ref( prm_pdb_atom );             // 13 - 16        Atom          name         Atom name.
	prm_os <<              prm_pdb_atom.get_alt_locn();                                    // 17             Character     altLoc       Alternate location indicator.
	prm_os <<              get_amino_acid_code_string( prm_pdb_atom );                     // 18 - 20        Residue name  resName      Residue name.
	prm_os << " ";
	prm_os << prm_res_id.get_chain_label();                                                // 22             Character     chainID      Chain identifier.
	                                                                                       // 23 - 26        Integer       resSeq       Residue sequence number.
	prm_os << setw( 5 ) << residue_name_with_insert_or_space;                              // 27             AChar         iCode        Code for insertion of residues.
	prm_os << "   ";
	prm_os << fixed     << setprecision( 3 );
	prm_os << setw( 8 ) << atom_coord.get_x();                                             // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
	prm_os << setw( 8 ) << atom_coord.get_y();                                             // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
	prm_os << setw( 8 ) << atom_coord.get_z();                                             // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
	prm_os << fixed     << setprecision( 2 );
	prm_os << setw( 6 ) << prm_pdb_atom.get_occupancy();                                   // 55 - 60        Real(6.2)     occupancy    Occupancy.
	prm_os << ( format( "%6.2f" ) % prm_pdb_atom.get_temp_factor() ).str().substr( 0, 6 ); // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
	const auto element_sym_strref = get_element_symbol_str_ref( prm_pdb_atom );            // 77 - 78        LString(2)    element      Element symbol, right-justified.
	const auto charge_strref      = get_charge_str_ref        ( prm_pdb_atom );            // 79 - 80        LString(2)    charge       Charge  on the atom.

	if ( ! element_sym_strref.empty() || ! charge_strref.empty() ) {
		prm_os << "          ";
		prm_os << right << setw( 2 ) << ( element_sym_strref.empty() ? "  "s : element_sym_strref.to_string() );
		if ( ! charge_strref.empty() ) {
			prm_os << charge_strref;
		}
	}

	// Output the result to the ostream and return it
	return prm_os;
}


/// \brief Generate a PDB file entry string for the specified residue_id and pdb_atom
///
/// \todo Change so that write_pdb_file_entry() calls this rather than vice versa
///
/// \relates pdb_atom
std::string cath::file::to_pdb_file_entry(const residue_id &prm_res_id,  ///< The residue_id to use in the PDB file entry string
                                          const pdb_atom   &prm_pdb_atom ///< The pdb_atom to use in the PDB file entry string
                                          ) {
	std::ostringstream output_ss;
	write_pdb_file_entry( output_ss, prm_res_id, prm_pdb_atom );
	return output_ss.str();
}


/// \brief TODOCUMENT
///
/// \relates pdb_atom
ostream & cath::file::operator<<(ostream        &prm_os,          ///< TODOCUMENT
                                 const pdb_atom &prm_pdb_atom ///< TODOCUMENT
                                 ) {
	prm_os << "Atom[";
	prm_os << prm_pdb_atom.get_element_type();
	prm_os << ", ";
	prm_os << prm_pdb_atom.get_coord();
	prm_os << "]";
	return prm_os;
}
