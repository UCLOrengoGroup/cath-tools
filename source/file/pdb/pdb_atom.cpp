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

#include "pdb_atom.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "file/pdb/pdb_base.h"
#include "structure/geometry/rotation.h"

#include <iomanip>
#include <sstream>

using namespace boost::math;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::algorithm::trim_copy;
using boost::algorithm::trim_left_copy;
using boost::algorithm::trim_right;
using boost::lexical_cast;
using boost::string_ref;

const string pdb_atom::PDB_ID_NITROGEN     ( "N"  );
const string pdb_atom::PDB_ID_CARBON_ALPHA ( "CA" );
const string pdb_atom::PDB_ID_CARBON       ( "C"  );

const string pdb_atom::PDB_ID_CARBON_BETA  ( "CB" );
const string pdb_atom::PDB_ID_OXYGEN       ( "O"  );

/// \brief TODOCUMENT
const pdb_record & pdb_atom::get_record_type() const {
	return record_type;
}

/// \brief TODOCUMENT
const size_t & pdb_atom::get_atom_serial() const {
	return atom_serial;
}

/// \brief TODOCUMENT
const string & pdb_atom::get_element_type_untrimmed() const {
	return the_element_type.get_element_type_untrimmed();
}

/// \brief TODOCUMENT
const string_ref & pdb_atom::get_element_type() const {
	return the_element_type.get_element_type();
}

/// \brief TODOCUMENT
const char & pdb_atom::get_alt_locn() const {
	return alt_locn;
}

/// \brief TODOCUMENT
const amino_acid & pdb_atom::get_amino_acid() const {
	return the_amino_acid;
}

/// \brief TODOCUMENT
const coord & pdb_atom::get_coord() const {
	return atom_coord;
}

/// \brief TODOCUMENT
const double & pdb_atom::get_occupancy() const {
	return occupancy;
}

/// \brief TODOCUMENT
const double & pdb_atom::get_temp_factor() const {
	return temp_factor;
}

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
char cath::file::get_amino_acid_letter(const pdb_atom &arg_pdb_atom ///< The pdb_atom to query
                                       ) {
	return arg_pdb_atom.get_amino_acid().get_letter();
}

/// \brief Get the three letter amino acid code (eg "SER") from a pdb_atom
string cath::file::get_amino_acid_code(const pdb_atom &arg_pdb_atom ///< The pdb_atom to query
                                       ) {
	return arg_pdb_atom.get_amino_acid().get_code();
}

/// \brief Get the three letter amino acid code (eg "SER") from a pdb_atom
string cath::file::get_amino_acid_name(const pdb_atom &arg_pdb_atom ///< The pdb_atom to query
                                       ) {
	return arg_pdb_atom.get_amino_acid().get_name();
}

/// \brief TODOCUMENT
///
/// \relates pdb_atom
ostream & cath::file::write_pdb_file_entry(ostream            &arg_os,          ///< TODOCUMENT
                                           const chain_label  &arg_chain_label, ///< TODOCUMENT
                                           const residue_name &arg_res_name,    ///< TODOCUMENT
                                           const pdb_atom     &arg_pdb_atom     ///< TODOCUMENT
                                           ) {
	// Sanity check the inputs
	if ( arg_res_name.get_is_null_residue_name() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Empty residue_name in cath::write_pdb_file_entry()"));
	}

	// Grab the necessary data
	const coord  atom_coord = arg_pdb_atom.get_coord();
	const string residue_name_with_insert_or_space = make_residue_name_string_with_insert_or_space( arg_res_name );

	// Output the entry to a temporary ostringstream
	// (to avoid needlessly messing around with the formatting options of the ostream)
	ostringstream atom_ss;
                                                                       // Comments with PDB format documentation
                                                                       // (http://www.wwpdb.org/documentation/format33/sect9.html#ATOM)
	atom_ss << left;
	atom_ss << setw( 6 ) << arg_pdb_atom.get_record_type();            //  1 -  6        Record name   "ATOM  " or "HETATM"
	atom_ss << right;
	atom_ss << setw( 5 ) << arg_pdb_atom.get_atom_serial();            //  7 - 11        Integer       serial       Atom  serial number.
	atom_ss << " ";
	atom_ss << setw( 4 ) << arg_pdb_atom.get_element_type_untrimmed(); // 13 - 16        Atom          name         Atom name.
	atom_ss <<              arg_pdb_atom.get_alt_locn();               // 17             Character     altLoc       Alternate location indicator.
	atom_ss <<              get_amino_acid_code( arg_pdb_atom );       // 18 - 20        Residue name  resName      Residue name.
	atom_ss << " ";
	atom_ss << arg_chain_label;                                        // 22             Character     chainID      Chain identifier.
	                                                                   // 23 - 26        Integer       resSeq       Residue sequence number.
	atom_ss << setw( 5 ) << residue_name_with_insert_or_space;         // 27             AChar         iCode        Code for insertion of residues.
	atom_ss << "   ";
	atom_ss << fixed     << setprecision( 3 );
	atom_ss << setw( 8 ) << atom_coord.get_x();                        // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
	atom_ss << setw( 8 ) << atom_coord.get_y();                        // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
	atom_ss << setw( 8 ) << atom_coord.get_z();                        // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
	atom_ss << fixed     << setprecision( 2 );
	atom_ss << setw( 6 ) << arg_pdb_atom.get_occupancy();              // 55 - 60        Real(6.2)     occupancy    Occupancy.
	atom_ss << setw( 6 ) << arg_pdb_atom.get_temp_factor();            // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
	                                                                   // 77 - 78        LString(2)    element      Element symbol, right-justified.
	                                                                   // 79 - 80        LString(2)    charge       Charge  on the atom.

	// Output the result to the ostream and return it
	arg_os << atom_ss.str();
	return arg_os;
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
