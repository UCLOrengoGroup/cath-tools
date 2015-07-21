/// \file
/// \brief The pdbs_acquirer class definitions

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

#include "pdbs_acquirer.h"

#include "common/clone/check_uptr_clone_against_this.h"
#include "common/logger.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"

using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace std;

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<pdbs_acquirer> pdbs_acquirer::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief TODOCUMENT
pdb_list_str_vec_pair pdbs_acquirer::get_pdbs_and_names(istream    &arg_istream,                ///< TODOCUMENT
                                                        const bool &arg_remove_partial_residues ///< TODOCUMENT
                                                        ) const {
	pair<pdb_list, str_vec> pdbs_and_names = do_get_pdbs_and_names( arg_istream );
	// Create a vector of PDBs to be superposed
	pdb_list &pdbs  = pdbs_and_names.first;
	str_vec  &names = pdbs_and_names.second;

	// Check the number of source files and then grab them
//	if (pdbs.size() < 2) {
//		logger::log_and_exit(
//			logger::return_code::TOO_FEW_PDBS_FOR_ALIGNMENT,
//			"ERROR: There aren't enough PDBs in the input to perform this alignment/superposition"
//		);
//	}

	// If the number of names doesn't match the number of PDBs then throw a wobbly
	if ( names.size() != pdbs.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The number of names doesn't match the number of PDBs"));
	}

	return arg_remove_partial_residues ? make_pair( pdb_list_of_backbone_complete_subset_pdbs( pdbs ), names )
	                                   : pdbs_and_names;
}
