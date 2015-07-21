/// \file
/// \brief The ostream_superposition_outputter class definitions

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

#include "ostream_superposition_outputter.h"

#include "common/clone/make_uptr_clone.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "superposition/superposition_context.h"
#include "superposition/superposition_io.h"

using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<superposition_outputter> ostream_superposition_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void ostream_superposition_outputter::do_output_superposition(const superposition_context &arg_superposition_context, ///< TODOCUMENT
		                                                      ostream                     &arg_ostream                ///< TODOCUMENT
		                                                      ) const {
	write_superposed_pdbs_to_ostream(
		arg_ostream,
		arg_superposition_context.get_superposition_cref(),
		arg_superposition_context.get_pdbs_cref(),
		false,
		true
	);
}

/// \brief TODOCUMENT
bool ostream_superposition_outputter::do_involves_display_spec() const {
	return false;
}

