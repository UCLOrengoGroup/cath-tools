/// \file
/// \brief The file_list_pdbs_acquirer class definitions

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

#include "file_list_pdbs_acquirer.h"

#include "common/clone/make_uptr_clone.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<pdbs_acquirer> file_list_pdbs_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
pdb_list_str_vec_pair file_list_pdbs_acquirer::do_get_pdbs_and_names(istream &/*arg_istream*/ ///< TODOCUMENT
                                                                     ) const {
	// Create a vector of PDBs to be superposed
	pdb_list pdbs;
	str_vec names;

	// Otherwise, load the PDBs from files
	const size_t num_input_files = files.size();
	for (size_t input_file_ctr = 0; input_file_ctr < num_input_files; ++input_file_ctr) {
		const path &input_filename = files[ input_file_ctr ];
		const pdb my_new_pdb = read_pdb_file( input_filename );
		pdbs.push_back( my_new_pdb );
		names.push_back( ( path( input_filename.stem() ) ).string() );
	}
	return make_pair( pdbs, names );
}

/// \brief Ctor for file_list_pdbs_acquirer.
file_list_pdbs_acquirer::file_list_pdbs_acquirer(const path_vec &arg_files ///< TODOCUMENT
                                                 ) : files( arg_files ) {
}


