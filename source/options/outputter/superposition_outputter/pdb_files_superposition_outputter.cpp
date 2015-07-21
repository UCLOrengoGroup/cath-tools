/// \file
/// \brief The pdb_files_superposition_outputter class definitions

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

#include "pdb_files_superposition_outputter.h"

#include "common/clone/make_uptr_clone.h"
#include "superposition/superposition_context.h"
#include "superposition/superposition_io.h"

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<superposition_outputter> pdb_files_superposition_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void pdb_files_superposition_outputter::do_output_superposition(const superposition_context &arg_superposition_context, ///< TODOCUMENT
                                                                ostream                     &/*arg_ostream*/            ///< TODOCUMENT
                                                                ) const {
	const superposition &the_superposition = arg_superposition_context.get_superposition_cref();
	const pdb_list      &pdbs              = arg_superposition_context.get_pdbs_cref();
	const str_vec       &names             = arg_superposition_context.get_names_cref();

	for (size_t pdb_ctr = 0; pdb_ctr < pdbs.size(); ++pdb_ctr) {
		const string &name                   = names[pdb_ctr];
		const path    superposition_output_file = output_dir / name;
		write_superposed_pdb_to_file(
			the_superposition,
			superposition_output_file.string(),
			pdbs[pdb_ctr],
			pdb_ctr
		);
	}
}

/// \brief TODOCUMENT
bool pdb_files_superposition_outputter::do_involves_display_spec() const {
	return false;
}

/// \brief Ctor for pdb_files_superposition_outputter.
pdb_files_superposition_outputter::pdb_files_superposition_outputter(const path &arg_output_dir
                                                                     ) : output_dir(arg_output_dir) {
}

