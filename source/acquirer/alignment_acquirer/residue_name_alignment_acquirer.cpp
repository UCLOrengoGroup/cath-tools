/// \file
/// \brief The residue_name_alignment_acquirer class definitions

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

#include "residue_name_alignment_acquirer.h"

#include <boost/numeric/conversion/cast.hpp>

#include "alignment/alignment.h"
#include "alignment/alignment_coord_extractor.h"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.h"
#include "alignment/residue_name_align/residue_name_aligner.h"
#include "alignment/residue_score/residue_scorer.h"
#include "common/clone/make_uptr_clone.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"
#include "structure/geometry/coord_list.h"
#include "structure/protein/protein.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"
#include "superposition/superpose_orderer.h"
#include "superposition/superposition.h"

#include <iostream>

using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

using boost::numeric_cast;

constexpr double residue_name_alignment_acquirer::RES_ALIGN_SCORE_CONSTANT;

/// \brief A standard do_clone method.
unique_ptr<alignment_acquirer> residue_name_alignment_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
pair<alignment, superpose_orderer> residue_name_alignment_acquirer::do_get_alignment_and_orderer(const pdb_list &arg_pdbs ///< TODOCUMENT
                                                                                                 ) const {
	const protein_list proteins_of_pdbs = build_protein_list_of_pdb_list( arg_pdbs );
	const size_t &num_pdbs = arg_pdbs.size();

	// Grab lists of the names of the residues in each PDB
	residue_name_vec_vec residue_names_of_pdbs;
	residue_names_of_pdbs.reserve( num_pdbs );
	for (const pdb &my_pdb : arg_pdbs) {
		residue_names_of_pdbs.push_back( my_pdb.get_residue_names_of_first_chain__backbone_unchecked() );

//		// *** TEMPORARY... ***
//		const str_vec residue_strings = my_pdb.get_residue_names_of_first_chain();
//		cerr << "Residue strings:";
//		for (const string &residue_string : residue_strings) {
//			cerr << " '" << residue_string << "'";
//		}
//		cerr << endl;
	}

	const alignment new_alignment = residue_name_align_and_residue_score_if_multi( residue_names_of_pdbs, residue_scorer(), proteins_of_pdbs );

	// For each pair, determine the number of residues in common and the RMSD
//	cerr << "Choosing spanning tree to use for superposing:" << endl;
	superpose_orderer my_orderer(num_pdbs);
	for (size_t pdb_ctr_1 = 0; pdb_ctr_1 < num_pdbs; ++pdb_ctr_1) {
		for (size_t pdb_ctr_2 = 0; pdb_ctr_2 < num_pdbs; ++pdb_ctr_2) {
			if (pdb_ctr_1 > pdb_ctr_2) {
				const pair<coord_list, coord_list> all_common_coords = alignment_coord_extractor::get_common_coords(
					new_alignment,
					arg_pdbs[ pdb_ctr_1 ],
					arg_pdbs[ pdb_ctr_2 ],
					common_residue_select_all_policy(),
					common_atom_select_ca_policy(),
					pdb_ctr_1,
					pdb_ctr_2
				);
				const size_t num_common_coords = all_common_coords.first.size();
				if (num_common_coords > MIN_NUM_COMMON_RESIDUES_TO_SUPERPOSE_PAIR) {
					const double rmsd = calc_pairwise_superposition_rmsd(all_common_coords.first, all_common_coords.second);
					const double score = numeric_cast<double>(num_common_coords) / (RES_ALIGN_SCORE_CONSTANT + rmsd);
//					cerr << "  - Comparing " << pdb_ctr_1;
//					cerr << "\tand " << pdb_ctr_2;
//					cerr << "\t--- num_common_coords : " << num_common_coords;
//					cerr << "\t, rmsd : " << rmsd;
//					cerr << "\t, score : " << score;
//					cerr << endl;
					my_orderer.set_score(pdb_ctr_1, pdb_ctr_2, score);
				}
			}
		}
	}

//	cerr << "Alignment built using residue names is : \n";
//	write_alignment_as_fasta_alignment( cerr, new_alignment, proteins_of_pdbs) << endl;
//	cerr << endl;

	// Return the results
	return make_pair(new_alignment, my_orderer);
}

