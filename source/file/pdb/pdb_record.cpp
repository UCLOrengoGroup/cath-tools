/// \file
/// \brief The pdb_record definitions

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

#include "pdb_record.h"

#include "exception/invalid_argument_exception.h"

#include <iostream>

using namespace cath::common;
using namespace cath::file;
using namespace std;

pdb_record cath::file::str_to_pdb_rec(const string &arg_string ///< TODOCUMENT
                                      ) {
	if      ( arg_string == "ATOM" || arg_string == "ATOM  " ) {
		return pdb_record::ATOM;
	}
	else if ( arg_string == "HETATM" ) {
		return pdb_record::HETATM;
	}
	else {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to recognise pdb_record type " + arg_string));
		return pdb_record::ATOM; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
	}
}

/// \brief Simple extraction operator for pdb_record
///
/// \relates pdb_record
istream & cath::file::operator>>(istream    &arg_is,        ///< The istream from which to extract the pdb_record
                                 pdb_record &arg_pdb_record ///< The pdb_record to populate
                                 ) {
	string input_string;
	arg_is >> input_string;
	arg_pdb_record = str_to_pdb_rec( input_string );
	return arg_is;
}

/// \brief TODOCUMENT
///
/// \relates pdb_record
ostream & cath::file::operator<<(ostream          &arg_os,        ///< TODOCUMENT
	                             const pdb_record &arg_pdb_record ///< TODOCUMENT
	                             ) {
	switch ( arg_pdb_record ) {
		case ( pdb_record::ATOM   ) : { arg_os << "ATOM";   break; }
		case ( pdb_record::HETATM ) : { arg_os << "HETATM"; break; }
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of pdb_record not recognised whilst inserting into an ostream"));
		}
	}
	return arg_os;
}
