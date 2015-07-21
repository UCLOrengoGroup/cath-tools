/// \file
/// \brief The ssap_ostream_alignment_outputter class definitions

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

#include "ssap_ostream_alignment_outputter.h"

#include "common/clone/make_uptr_clone.h"
#include "alignment/alignment_context.h"
#include "alignment/alignment_io.h"
#include "alignment/pair_alignment.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_list.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

using namespace cath::align;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<alignment_outputter> ssap_ostream_alignment_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void ssap_ostream_alignment_outputter::do_output_alignment(const alignment_context &arg_alignment_context, ///< TODOCUMENT
                                                           ostream                 &arg_ostream            ///< TODOCUMENT
                                                           ) const {
	const alignment the_alignment = arg_alignment_context.get_alignment();
	check_alignment_is_a_pair( the_alignment );

	const protein_list temp_protein_list = build_protein_list_of_pdb_list_and_names(
		arg_alignment_context.get_pdbs(),
		arg_alignment_context.get_names()
	);
	output_alignment_to_cath_ssap_legacy_format(
		arg_ostream,
		the_alignment,
		temp_protein_list[ 0 ],
		temp_protein_list[ 1 ]
	);
}

/// \brief TODOCUMENT
bool ssap_ostream_alignment_outputter::do_involves_display_spec() const {
	return false;
}

