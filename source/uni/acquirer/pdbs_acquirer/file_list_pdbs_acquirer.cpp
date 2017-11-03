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

#include "file_list_pdbs_acquirer.hpp"

#include "common/boost_addenda/range/indices.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "file/name_set/name_set_list.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"

using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;

using boost::filesystem::path;
using std::istream;
using std::make_pair;
using std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<pdbs_acquirer> file_list_pdbs_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
pdb_list_name_set_list_pair file_list_pdbs_acquirer::do_get_pdbs_and_names(istream &/*arg_istream*/ ///< TODOCUMENT
                                                                           ) const {
	// Create a vector of PDBs to be superposed
	pdb_list pdbs;
	name_set_vec names;

	// Otherwise, load the PDBs from files
	const size_t num_input_files = files.size();
	for (const size_t &input_file_ctr : indices( num_input_files ) ) {
		const path &input_filename = files[ input_file_ctr ];
		const pdb my_new_pdb = read_pdb_file( input_filename );
		pdbs.push_back( my_new_pdb );
		names.emplace_back( input_filename );
	}
	return make_pair( pdbs, name_set_list{ std::move( names ) } );
}

/// \brief Ctor for file_list_pdbs_acquirer.
file_list_pdbs_acquirer::file_list_pdbs_acquirer(path_vec arg_files ///< TODOCUMENT
                                                 ) : files { std::move( arg_files ) } {
}


