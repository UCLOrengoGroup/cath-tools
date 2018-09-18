/// \file
/// \brief The ssap_ostream_alignment_outputter class definitions

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

#include "ssap_ostream_alignment_outputter.hpp"

#include "alignment/alignment_context.hpp"
#include "alignment/io/alignment_io.hpp"
#include "alignment/pair_alignment.hpp"
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
using namespace cath::chop;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<alignment_outputter> ssap_ostream_alignment_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void ssap_ostream_alignment_outputter::do_output_alignment(const alignment_context &prm_alignment_context, ///< TODOCUMENT
                                                           ostream                 &prm_ostream            ///< TODOCUMENT
                                                           ) const {
	const region_vec_opt_vec &regions_opts  = get_regions( prm_alignment_context );
	const alignment          &the_alignment = prm_alignment_context.get_alignment();
	check_alignment_is_a_pair( the_alignment );

	const protein_list temp_protein_list = build_protein_list_of_pdb_list_and_names(
		get_pdbs      ( prm_alignment_context ),
		get_name_sets ( prm_alignment_context )
	);

	output_alignment_to_cath_ssap_legacy_format(
		prm_ostream,
		the_alignment,
		temp_protein_list[ 0 ],
		temp_protein_list[ 1 ],
		regions_opts     [ 0 ],
		regions_opts     [ 1 ]
	);
}

/// \brief TODOCUMENT
bool ssap_ostream_alignment_outputter::do_involves_display_spec() const {
	return false;
}

/// \brief Get a name for this alignment_outputter
string ssap_ostream_alignment_outputter::do_get_name() const {
	return "SSAP";
}
