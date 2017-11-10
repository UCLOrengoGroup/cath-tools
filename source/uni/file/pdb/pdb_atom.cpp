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

#include "structure/geometry/rotation.hpp"

#include <iomanip>
#include <sstream>

using namespace boost::math;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::format;
using boost::io::ios_flags_saver;

constexpr size_t pdb_atom::MIN_NUM_PDB_COLS;
constexpr size_t pdb_atom::MAX_NUM_PDB_COLS;

/// \brief TODOCUMENT
void pdb_atom::rotate(const rotation &arg_rotation ///< TODOCUMENT
                      ) {
	cath::geom::rotate( arg_rotation, atom_coord );
}

/// \brief TODOCUMENT
void pdb_atom::operator+=(const coord &arg_coord ///< TODOCUMENT
                          ) {
	atom_coord += arg_coord;
}

/// \brief TODOCUMENT
void pdb_atom::operator-=(const coord &arg_coord ///< TODOCUMENT
                          ) {
	atom_coord -= arg_coord;
}

/// \brief Get the one letter amino acid code (eg 'S') from a pdb_atom
char cath::file::get_amino_acid_letter_tolerantly(const pdb_atom &arg_pdb_atom ///< The pdb_atom to query
                                                  ) {
	return arg_pdb_atom.get_amino_acid().get_letter_tolerantly();
}

/// \brief Get the three letter amino acid code (eg "SER") char_3_arr from a pdb_atom
char_3_arr cath::file::get_amino_acid_code(const pdb_atom &arg_pdb_atom ///< The pdb_atom to query
                                           ) {
	return arg_pdb_atom.get_amino_acid().get_code();
}

/// \brief Get the three letter amino acid code (eg "SER") string from a pdb_atom
string cath::file::get_amino_acid_code_string(const pdb_atom &arg_pdb_atom ///< The pdb_atom to query
                                       ) {
	return get_code_string( arg_pdb_atom.get_amino_acid() );
}

/// \brief Get the three letter amino acid code (eg "SER") from a pdb_atom
string cath::file::get_amino_acid_name(const pdb_atom &arg_pdb_atom ///< The pdb_atom to query
                                       ) {
	return arg_pdb_atom.get_amino_acid().get_name();
}

/// \brief TODOCUMENT
///
/// \relates pdb_atom
ostream & cath::file::write_pdb_file_entry(ostream          &arg_os,      ///< TODOCUMENT
                                           const residue_id &arg_res_id,  ///< TODOCUMENT
                                           const pdb_atom   &arg_pdb_atom ///< TODOCUMENT
                                           ) {
	// Sanity check the inputs
	if ( arg_res_id.get_residue_name().is_null() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Empty residue_name in cath::write_pdb_file_entry()"));
	}

	// Save the state of the ostream's flags (and reset them in this guard's dtor)
	const ios_flags_saver stream_flags_guard{ arg_os };

	// Grab the necessary data
	const coord  atom_coord                        = arg_pdb_atom.get_coord();
	const string residue_name_with_insert_or_space = make_residue_name_string_with_insert_or_space(
		arg_res_id.get_residue_name()
	);

	                                                                                       // Comments with PDB format documentation
	                                                                                       // (http://www.wwpdb.org/documentation/format33/sect9.html#ATOM)
	arg_os << left;
	arg_os << setw( 6 ) << arg_pdb_atom.get_record_type();                                 //  1 -  6        Record name   "ATOM  " or "HETATM"
	arg_os << right;
	arg_os << setw( 5 ) << arg_pdb_atom.get_atom_serial();                                 //  7 - 11        Integer       serial       Atom  serial number.
	arg_os << " ";
	arg_os << setw( 4 ) << get_element_type_untrimmed_str_ref( arg_pdb_atom );             // 13 - 16        Atom          name         Atom name.
	arg_os <<              arg_pdb_atom.get_alt_locn();                                    // 17             Character     altLoc       Alternate location indicator.
	arg_os <<              get_amino_acid_code_string( arg_pdb_atom );                     // 18 - 20        Residue name  resName      Residue name.
	arg_os << " ";
	arg_os << arg_res_id.get_chain_label();                                                // 22             Character     chainID      Chain identifier.
	                                                                                       // 23 - 26        Integer       resSeq       Residue sequence number.
	arg_os << setw( 5 ) << residue_name_with_insert_or_space;                              // 27             AChar         iCode        Code for insertion of residues.
	arg_os << "   ";
	arg_os << fixed     << setprecision( 3 );
	arg_os << setw( 8 ) << atom_coord.get_x();                                             // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
	arg_os << setw( 8 ) << atom_coord.get_y();                                             // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
	arg_os << setw( 8 ) << atom_coord.get_z();                                             // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
	arg_os << fixed     << setprecision( 2 );
	arg_os << setw( 6 ) << arg_pdb_atom.get_occupancy();                                   // 55 - 60        Real(6.2)     occupancy    Occupancy.
	arg_os << ( format( "%6.2f" ) % arg_pdb_atom.get_temp_factor() ).str().substr( 0, 6 ); // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
	const auto element_sym_strref = get_element_symbol_str_ref( arg_pdb_atom );            // 77 - 78        LString(2)    element      Element symbol, right-justified.
	const auto charge_strref      = get_charge_str_ref        ( arg_pdb_atom );            // 79 - 80        LString(2)    charge       Charge  on the atom.

	if ( ! element_sym_strref.empty() || ! charge_strref.empty() ) {
		arg_os << "          ";
		arg_os << right << setw( 2 ) << ( element_sym_strref.empty() ? "  "s : element_sym_strref.to_string() );
		if ( ! charge_strref.empty() ) {
			arg_os << charge_strref;
		}
	}

	// Output the result to the ostream and return it
	return arg_os;
}


/// \brief Generate a PDB file entry string for the specified residue_id and pdb_atom
///
/// \todo Change so that write_pdb_file_entry() calls this rather than vice versa
///
/// \relates pdb_atom
std::string cath::file::to_pdb_file_entry(const residue_id &arg_res_id,  ///< The residue_id to use in the PDB file entry string
                                          const pdb_atom   &arg_pdb_atom ///< The pdb_atom to use in the PDB file entry string
                                          ) {
	std::ostringstream output_ss;
	write_pdb_file_entry( output_ss, arg_res_id, arg_pdb_atom );
	return output_ss.str();
}


/// \brief TODOCUMENT
///
/// \relates pdb_atom
ostream & cath::file::operator<<(ostream        &arg_os,          ///< TODOCUMENT
                                 const pdb_atom &arg_pdb_atom ///< TODOCUMENT
                                 ) {
	arg_os << "Atom[";
	arg_os << arg_pdb_atom.get_element_type();
	arg_os << ", ";
	arg_os << arg_pdb_atom.get_coord();
	arg_os << "]";
	return arg_os;
}
