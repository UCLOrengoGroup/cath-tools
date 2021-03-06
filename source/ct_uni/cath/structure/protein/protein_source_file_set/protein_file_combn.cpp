/// \file
/// \brief The protein_file_combn definitions

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

#include "protein_file_combn.hpp"

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_pdb.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_pdb_and_calc.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_pdb_and_dssp_and_calc.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_pdb_dssp_and_sec.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_from_wolf_and_sec.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::std;

/// \brief Convert a protein_file_combn to (a smart pointer to) the equivalent protein_source_file_set
///
/// \relates protein_file_combn
unique_ptr<const protein_source_file_set> cath::get_protein_source_file_set(const protein_file_combn &prm_protein_file_combn /// The protein_file_combn to be converted
                                                                            ) {
	switch ( prm_protein_file_combn ) {
		case( protein_file_combn::WOLF_SEC          ) : { return { ::std::make_unique< protein_from_wolf_and_sec          >() }; break; }
		case( protein_file_combn::PDB               ) : { return { ::std::make_unique< protein_from_pdb                   >() }; break; }
		case( protein_file_combn::PDB_DSSP_SEC      ) : { return { ::std::make_unique< protein_from_pdb_dssp_and_sec      >() }; break; }
		case( protein_file_combn::PDB_DSSP_AND_CALC ) : { return { ::std::make_unique< protein_from_pdb_and_dssp_and_calc >() }; break; }
		case( protein_file_combn::PDB_AND_CALC      ) : { return { ::std::make_unique< protein_from_pdb_and_calc          >() }; break; }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("protein_file_combn is not recognised"));
}

/// \brief Simple extraction operator for protein_file_combn
///
/// \relates protein_file_combn
istream & cath::operator>>(istream            &prm_is,                ///< The istream from which to extract the protein_file_combn
                           protein_file_combn &prm_protein_file_combn ///< The protein_file_combn to populate
                           ) {
	string input_string;
	prm_is >> input_string;
	if ( input_string == "WOLF_SEC" ) {
		prm_protein_file_combn = protein_file_combn::WOLF_SEC;
	}
	else if ( input_string == "PDB_SIMPLE" ) {
		prm_protein_file_combn = protein_file_combn::PDB;
	}
	else if ( input_string == "PDB_DSSP_SEC" ) {
		prm_protein_file_combn = protein_file_combn::PDB_DSSP_SEC;
	}
	else if ( input_string == "PDB_DSSP" ) {
		prm_protein_file_combn = protein_file_combn::PDB_DSSP_AND_CALC;
	}
	else if ( input_string == "PDB" ) {
		prm_protein_file_combn = protein_file_combn::PDB_AND_CALC;
	}
	else {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to recognise protein_file_combn type " + input_string));
	}
	return prm_is;
}

/// \brief Simple insertion operator for protein_file_combn
///
/// \relates protein_file_combn
ostream & cath::operator<<(ostream                  &prm_os,                ///< The stream to which to output
                           const protein_file_combn &prm_protein_file_combn ///< The protein_file_combn to output to the stream
                           ) {
	switch( prm_protein_file_combn ) {
		case( protein_file_combn::WOLF_SEC ) : {
			prm_os << "WOLF_SEC";
			break;
		}
		case( protein_file_combn::PDB ) : {
			prm_os << "PDB_SIMPLE";
			break;
		}
		case( protein_file_combn::PDB_DSSP_SEC ) : {
			prm_os << "PDB_DSSP_SEC";
			break;
		}
		case( protein_file_combn::PDB_DSSP_AND_CALC ) : {
			prm_os << "PDB_DSSP";
			break;
		}
		case( protein_file_combn::PDB_AND_CALC ) : {
			prm_os << "PDB";
			break;
		}
	}
	return prm_os;
}

