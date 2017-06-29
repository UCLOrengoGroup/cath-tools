/// \file
/// \brief The fasta_ostream_alignment_outputter class definitions

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

#include "fasta_ostream_alignment_outputter.hpp"

#include "alignment/alignment_context.hpp"
#include "alignment/io/alignment_io.hpp"
#include "chopping/region/region.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath::align;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<alignment_outputter> fasta_ostream_alignment_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void fasta_ostream_alignment_outputter::do_output_alignment(const alignment_context &arg_alignment_context, ///< TODOCUMENT
                                                            ostream                 &arg_ostream            ///< TODOCUMENT
                                                            ) const {
	write_alignment_as_fasta_alignment(
		arg_ostream,
		arg_alignment_context.get_alignment(),
		build_protein_list_of_pdb_list_and_names(
			get_restricted_pdbs( arg_alignment_context ),
			get_name_sets      ( arg_alignment_context )
		)
	);
}

/// \brief TODOCUMENT
bool fasta_ostream_alignment_outputter::do_involves_display_spec() const {
	return false;
}

