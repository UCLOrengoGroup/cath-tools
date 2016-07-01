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

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "common/algorithm/contains.h"
#include "common/string/sub_string_parser.h"
#include "exception/invalid_argument_exception.h"
#include "structure/geometry/rotation.h"
#include "file/pdb/pdb_base.h"

#include <iomanip>
#include <iostream> // *** TEMPORARY ***
#include <sstream>

using namespace boost::math;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::algorithm::starts_with;
using boost::algorithm::trim_copy;
using boost::algorithm::trim_left_copy;
using boost::algorithm::trim_right;
using boost::lexical_cast;

const string pdb_atom::PDB_ID_NITROGEN     ( "N"  );
const string pdb_atom::PDB_ID_CARBON_ALPHA ( "CA" );
const string pdb_atom::PDB_ID_CARBON       ( "C"  );

const string pdb_atom::PDB_ID_CARBON_BETA  ( "CB" );
const string pdb_atom::PDB_ID_OXYGEN       ( "O"  );

/// \brief Ctor for pdb_atom
pdb_atom::pdb_atom(const pdb_record &arg_record_type,  ///< TODOCUMENT
                   const size_t     &arg_atom_serial,  ///< TODOCUMENT
                   const string     &arg_element_type, ///< TODOCUMENT
                   const char       &arg_alt_locn,     ///< TODOCUMENT
                   const amino_acid &arg_amino_acid,   ///< TODOCUMENT
                   const coord      &arg_coord,        ///< TODOCUMENT
                   const double     &arg_occupancy,    ///< TODOCUMENT
                   const double     &arg_temp_factor   ///< TODOCUMENT
                   ) : record_type           ( arg_record_type               ),
                       atom_serial           ( arg_atom_serial               ),
                       element_type_untrimmed( arg_element_type              ),
                       element_type          ( trim_copy( arg_element_type ) ),
                       alt_locn              ( arg_alt_locn                  ),
                       the_amino_acid        ( arg_amino_acid                ),
                       atom_coord            ( arg_coord                     ),
                       occupancy             ( arg_occupancy                 ),
                       temp_factor           ( arg_temp_factor               ) {
	using boost::math::isfinite;

	if ( ! isfinite( occupancy ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Argument occupancy must be a normal, finite floating-point number"));
	}
	if ( ! isfinite( temp_factor ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Argument temp_factor must be a normal, finite floating-point number"));
	}
}

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
	return element_type_untrimmed;
}

/// \brief TODOCUMENT
const string & pdb_atom::get_element_type() const {
	return element_type;
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

/// \brief Return whether this line purports to be an record of the specified pdb_record
///
/// Note: This does NOT check whether the line is a valid record;
///       Use atom_record_parse_problem() for that.
bool cath::file::is_pdb_record_of_type(const string     &arg_pdb_record_string, ///< The string to check
                                       const pdb_record &arg_record_type        ///< TODOCUMENT
                                       ) {
	if ( arg_pdb_record_string.empty() ) {
		return false;
	}
	switch ( arg_record_type ) {
		case ( pdb_record::ATOM   ) : { return starts_with( arg_pdb_record_string, "ATOM  " ); }
		case ( pdb_record::HETATM ) : { return starts_with( arg_pdb_record_string, "HETATM" ); }
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of pdb_record not recognised whilst checking is_pdb_record_of_type()"));
			return false; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
}

/// \brief Convert the specified three-letter string and pdb_record to an amino_acid
///
/// This is less strict than the amino_acid ctor because it allows certain combinations
/// to be accepted for decay to UNK/X. See function body for the accepted combinations.
///
/// \relates pdb_atom
amino_acid cath::file::get_amino_acid_of_string_and_record(const string     &arg_aa_string,  ///< The three-letter amino acid string (eg "SER")
                                                           const pdb_record &arg_record_type ///< Whether this is an ATOM or HETATM record
                                                           ) {
	if ( arg_aa_string.length() != 3 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Amino acid string must contain three characters"));
	}
	try {
		return { arg_aa_string, arg_record_type };
	}
	catch (...) {
		const set<pair<string, pdb_record>> pdb_aa_entries_allowed_to_decay_to_unk = {
			{ "MLY", pdb_record::HETATM },
			{ "MSE", pdb_record::ATOM   },
			{ "MSE", pdb_record::HETATM },
		};
		if ( contains( pdb_aa_entries_allowed_to_decay_to_unk, make_pair( arg_aa_string, arg_record_type ) ) ) {
			return { "UNK", arg_record_type };
		}
		throw;
	}
}

/// \brief Parse the amino_acid from the specified PDB atom record string using the
///        looser pdb_atom criteria that allow a few extra amino acids, which decay to UNK/X.
///
/// See the body of get_amino_acid_of_string_and_record() for the accepted combinations.
///
/// \relates pdb_atom
amino_acid cath::file::parse_amino_acid_from_pdb_atom_record(const string &arg_pdb_atom_record_string ///< TODOCUMENT
                                                             ) {
	sub_string_parser p;
	const pdb_record       record_type = str_to_pdb_rec( p.substr_as_str_ref( arg_pdb_atom_record_string,         0, 6  )); //  1 -  6        Record name   "ATOM  "
	const string           the_a_a     =                                      arg_pdb_atom_record_string.substr( 17, 3  ) ; // 18 - 20        Residue name  resName      Residue name.
	return get_amino_acid_of_string_and_record(
		the_a_a,
		record_type
	);
}

/// \brief Return a string containing the parse problem with a PDB ATOM/HETATM record string or "" if no problem
///
/// \relates pdb_atom
///
/// Note: This does NOT check whether the line is a valid record.
///       Use atom_record_parse_problem() for that.
std::pair<pdb_atom_parse_status, std::string> cath::file::pdb_record_parse_problem(const string &arg_pdb_atom_record_string ///< The string to check
                                                                                   ) {
	if ( ! is_pdb_record_of_type( arg_pdb_atom_record_string, pdb_record::ATOM ) && ! is_pdb_record_of_type( arg_pdb_atom_record_string, pdb_record::HETATM ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot check for atom record parse problems because string is not an ATOM record"));
	}
	if ( arg_pdb_atom_record_string.length() > pdb_base::MAX_NUM_PDB_COLS ) {
		return { pdb_atom_parse_status::ABORT, "Is too long" };
	}
	if ( arg_pdb_atom_record_string.at( 11 ) != ' '   ) {
		return { pdb_atom_parse_status::ABORT, "Does not contain a space at column 12" };
	}
	if ( arg_pdb_atom_record_string.at( 20 ) != ' '   ) {
		return { pdb_atom_parse_status::ABORT, "Does not contain a space at column 21" };
	}
	if ( arg_pdb_atom_record_string.at( 27 ) != ' ' || arg_pdb_atom_record_string.at( 28 ) != ' ' || arg_pdb_atom_record_string.at( 29 ) != ' ' ) {
		return { pdb_atom_parse_status::ABORT, "Does not contain spaces at columns 28-30" };
	}
//	if ( arg_pdb_atom_record_string.at( 16 ) != ' ' && arg_pdb_atom_record_string.at( 16 ) != 'A' ) {
//		return { pdb_atom_parse_status::SKIP, "Has alternate location indicator other than 'A' or' '" };
//	}
	try {
		parse_amino_acid_from_pdb_atom_record( arg_pdb_atom_record_string );
	}
	catch (...) {
		return { pdb_atom_parse_status::SKIP, "Do not recognise amino acid entry: " + arg_pdb_atom_record_string.substr( 17, 3 ) };
	}
	return { pdb_atom_parse_status::OK, "" };
}

/// \brief TODOCUMENT
///
/// \relates pdb_atom
chain_resname_atom_tuple cath::file::parse_pdb_atom_record(string arg_pdb_atom_record_string ///< TODOCUMENT
                                                           ) {
	if ( pdb_record_parse_problem( arg_pdb_atom_record_string ).first != pdb_atom_parse_status::OK ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse string - is not a valid PDB ATOM/HETATM record"));
	}

	const amino_acid the_amino_acid = parse_amino_acid_from_pdb_atom_record( arg_pdb_atom_record_string );

	sub_string_parser p;

	if ( isspace( arg_pdb_atom_record_string.at( 25 ) ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("PDB ATOM/HETATOM record malformed: space found in column 26 which should contain the end of a right justified residue number"));
	}

	try {
		                                                                                                                        // Comments with PDB format documentation
		                                                                                                                        // (http://www.wwpdb.org/documentation/format33/sect9.html#ATOM)
		const pdb_record       record_type = str_to_pdb_rec( p.substr_as_str_ref( arg_pdb_atom_record_string,         0, 6  )); //  1 -  6        Record name   "ATOM  "
		const size_t           serial      =                  p.substr_as_size_t( arg_pdb_atom_record_string,         6, 5  ) ; //  7 - 11        Integer       serial       Atom  serial number.
		const string           element     =                                      arg_pdb_atom_record_string.substr( 12, 4  ) ; // 13 - 16        Atom          name         Atom name.
		const char            &alt_locn    =                                      arg_pdb_atom_record_string.at(     16     ) ; // 17             Character     altLoc       Alternate location indicator.
//		const string           the_a_a     =                                      arg_pdb_atom_record_string.substr( 17, 3  ) ; // 18 - 20        Residue name  resName      Residue name.
		const char            &chain_char  =                                      arg_pdb_atom_record_string.at(     21     ) ; // 22             Character     chainID      Chain identifier.
		const int              res_num     =                     p.substr_as_int( arg_pdb_atom_record_string,        22, 4  ) ; // 23 - 26        Integer       resSeq       Residue sequence number.
		const char            &insert_code =                                      arg_pdb_atom_record_string.at    ( 26     ) ; // 27             AChar         iCode        Code for insertion of residues.
		const double           coord_x     =                  p.substr_as_double( arg_pdb_atom_record_string,        30, 8  ) ; // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		const double           coord_y     =                  p.substr_as_double( arg_pdb_atom_record_string,        38, 8  ) ; // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		const double           coord_z     =                  p.substr_as_double( arg_pdb_atom_record_string,        46, 8  ) ; // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		const double           occupancy   =                  p.substr_as_double( arg_pdb_atom_record_string,        54, 6  ) ; // 55 - 60        Real(6.2)     occupancy    Occupancy.
		const double           temp_factor =                  p.substr_as_double( arg_pdb_atom_record_string,        60, 6  ) ; // 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
//		const string           element     =                           trim_copy( arg_pdb_atom_record_string.substr( 76, 3 ));  // 77 - 78        LString(2)    element      Element symbol, right-justified.
//		const string           charge      =                           trim_copy( arg_pdb_atom_record_string.substr( 78, 2 ));  // 79 - 80        LString(2)    charge       Charge  on the atom.

		return make_tuple(
			chain_label( chain_char ),
			make_residue_name_with_non_insert_char( res_num, insert_code, ' ' ),
			pdb_atom(
				record_type,
				serial,
				element,
				alt_locn,
				the_amino_acid,
				coord( coord_x, coord_y, coord_z ),
				occupancy,
				temp_factor
			)
		);
	}
	catch (const boost::bad_lexical_cast &) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to cast a column whilst parsing a PDB ATOM record, which probably means it's malformed.\nRecord was \"" + arg_pdb_atom_record_string + "\""));
	}
	catch (const invalid_argument &) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to cast a column whilst parsing a PDB ATOM record, which probably means it's malformed.\nRecord was \"" + arg_pdb_atom_record_string + "\""));
	}
	catch (const out_of_range &) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Casted column out of range whilst parsing a PDB ATOM record, which probably means it's malformed.\nRecord was \"" + arg_pdb_atom_record_string + "\""));
	}
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
