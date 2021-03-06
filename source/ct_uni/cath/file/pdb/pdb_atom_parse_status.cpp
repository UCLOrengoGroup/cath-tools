/// \file
/// \brief The pdb_atom_parse_status definitions

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

#include "pdb_atom_parse_status.hpp"

#include "cath/common/exception/invalid_argument_exception.hpp"

#include <iostream>

using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::std;

//pdb_atom_parse_status cath::file::str_to_pdb_rec(const string &prm_string ///< TODOCUMENT
//                                                 ) {
//	if      ( prm_string == "OK" ) {
//		return pdb_atom_parse_status::OK;
//	}
//	else if ( prm_string == "SKIP" ) {
//			return pdb_atom_parse_status::SKIP;
//		}
//	else if ( prm_string == "ABORT" ) {
//		return pdb_atom_parse_status::ABORT;
//	}
//	else {
//		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to recognise pdb_atom_parse_status type " + prm_string));
//		return pdb_atom_parse_status::ABORT; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
//	}
//}

///// \brief Simple extraction operator for pdb_atom_parse_status
/////
///// \relates pdb_atom_parse_status
//istream & cath::file::operator>>(istream               &prm_is,                   ///< The istream from which to extract the pdb_atom_parse_status
//                                 pdb_atom_parse_status &prm_pdb_atom_parse_status ///< The pdb_atom_parse_status to populate
//                                 ) {
//	string input_string;
//	prm_is >> input_string;
//	prm_pdb_atom_parse_status = str_to_pdb_rec( input_string );
//	return prm_is;
//}

///// \brief TODOCUMENT
/////
///// \relates pdb_atom_parse_status
//ostream & cath::file::operator<<(ostream                     &prm_os,                   ///< TODOCUMENT
//	                               const pdb_atom_parse_status &prm_pdb_atom_parse_status ///< TODOCUMENT
//	                               ) {
//	switch ( prm_pdb_atom_parse_status ) {
//		case ( pdb_atom_parse_status::OK    ) : { prm_os << "OK"   ; return prm_os ; }
//		case ( pdb_atom_parse_status::SKIP  ) : { prm_os << "SKIP" ; return prm_os ; }
//		case ( pdb_atom_parse_status::ABORT ) : { prm_os << "ABORT"; return prm_os ; }
//	}
//	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of pdb_atom_parse_status not recognised whilst inserting into an ostream"));
//}
