/// \file
/// \brief The istream_pdbs_acquirer class definitions

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

#include "istream_pdbs_acquirer.hpp"

#include <boost/lexical_cast.hpp>

#include "common/clone/make_uptr_clone.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"

using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace std;

using boost::lexical_cast;

/// \brief A standard do_clone method.
unique_ptr<pdbs_acquirer> istream_pdbs_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
pdb_list_str_vec_pair istream_pdbs_acquirer::do_get_pdbs_and_names(istream &arg_istream ///< TODOCUMENT
                                                                   ) const {
	// Create a vector of PDBs to be superposed
	pdb_list pdbs;
	str_vec names;

	// Read PDBs from the_istream
	pdbs = read_end_separated_pdb_files( arg_istream );
	names.assign( pdbs.size(), string() );
	for (size_t names_ctr = 0; names_ctr < names.size(); ++names_ctr) {
		names[ names_ctr ] = "PDB_"
		                     + lexical_cast<string>( names_ctr + 1 )
		                     + "_from_stdin (with "
		                     + lexical_cast<string>( pdbs[ names_ctr ].get_num_atoms() )
		                     + " atoms)";
	}
	return make_pair( pdbs, names );
}

