/// \file
/// \brief The pdb_record definitions

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

#include "pdb_record.hpp"

#include "exception/invalid_argument_exception.hpp"

#include <iostream>

using namespace cath::common;
using namespace cath::file;
using namespace std;

/// \brief Simple extraction operator for pdb_record
///
/// \relates pdb_record
istream & cath::file::operator>>(istream    &arg_is,        ///< The istream from which to extract the pdb_record
                                 pdb_record &arg_pdb_record ///< The pdb_record to populate
                                 ) {
	string input_string;
	arg_is >> input_string;
	arg_pdb_record = pdb_rec_of_str( input_string );
	return arg_is;
}

/// \brief TODOCUMENT
///
/// \relates pdb_record
ostream & cath::file::operator<<(ostream          &arg_os,        ///< TODOCUMENT
	                             const pdb_record &arg_pdb_record ///< TODOCUMENT
	                             ) {
	switch ( arg_pdb_record ) {
		case ( pdb_record::ATOM   ) : { arg_os << "ATOM"   ; return arg_os ; }
		case ( pdb_record::HETATM ) : { arg_os << "HETATM" ; return arg_os ; }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of pdb_record not recognised whilst inserting into an ostream"));
}
