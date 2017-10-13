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

#include "common/algorithm/transform_build.hpp"
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
using namespace std;

using boost::lexical_cast;

/// \brief A standard do_clone method.
unique_ptr<pdbs_acquirer> istream_pdbs_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
pdb_list_name_set_list_pair istream_pdbs_acquirer::do_get_pdbs_and_names(istream &arg_istream ///< TODOCUMENT
                                                                         ) const {
	using std::to_string;

	// Read PDBs from the_istream
	const pdb_list pdbs = read_end_separated_pdb_files( arg_istream );

	return make_pair(
		pdbs,
		name_set_list{
			transform_build<name_set_vec>(
				indices( pdbs.size() ),
				[&] (const size_t &x) {
					return name_set{
						"PDB_"
						+ to_string( x + 1 )
						+ "_from_stdin (with "
						+ to_string( pdbs[ x ].get_num_atoms() )
						+ " atoms)"
					};
				}
			)
		}
	);
}

