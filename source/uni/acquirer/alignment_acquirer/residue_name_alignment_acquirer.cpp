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

#include "residue_name_alignment_acquirer.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include "alignment/alignment.hpp"
#include "alignment/alignment_coord_extractor.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "alignment/residue_name_align/residue_name_aligner.hpp"
#include "alignment/residue_score/residue_scorer.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/graph/spanning_tree.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "file/strucs_context.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "superposition/superposition.hpp"

#include <iostream>

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
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

/// \brief Specify that this does require backbone-complete input
bool residue_name_alignment_acquirer::do_requires_backbone_complete_input() const {
	return true;
}

/// \brief TODOCUMENT
pair<alignment, size_size_pair_vec> residue_name_alignment_acquirer::do_get_alignment_and_spanning_tree(const strucs_context &prm_strucs_context ///< TODOCUMENT
                                                                                                        ) const {
	const protein_list  proteins_of_pdbs = build_protein_list( prm_strucs_context );
	const auto         &the_pdbs         = prm_strucs_context.get_pdbs();
	const size_t        num_pdbs         = the_pdbs.size();

	// Grab lists of the names of the residues in each PDB
	residue_name_vec_vec residue_names_of_pdbs;
	residue_names_of_pdbs.reserve( num_pdbs );
	for (const pdb &my_pdb : the_pdbs) {
		residue_names_of_pdbs.push_back(
			transform_build<residue_name_vec>(
				my_pdb.get_residue_ids_of_first_chain__backbone_unchecked(),
				[] (const residue_id &x) { return x.get_residue_name(); }
			)
		);

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
	size_size_doub_tpl_vec edges;
	for (const size_t &pdb_ctr_1 : indices( num_pdbs ) ) {
		for (const size_t &pdb_ctr_2 : indices( pdb_ctr_1 ) ) {
			const pair<coord_list, coord_list> all_common_coords = alignment_coord_extractor::get_common_coords(
				new_alignment,
				the_pdbs[ pdb_ctr_1 ],
				the_pdbs[ pdb_ctr_2 ],
				common_residue_select_all_policy(),
				common_atom_select_ca_policy(),
				pdb_ctr_1,
				pdb_ctr_2
			);
			const size_t num_common_coords = all_common_coords.first.size();
			if ( num_common_coords > MIN_NUM_COMMON_RESIDUES_TO_SUPERPOSE_PAIR ) {
				const double rmsd  = calc_pairwise_superposition_rmsd( all_common_coords.first, all_common_coords.second );
				const double score = numeric_cast<double>( num_common_coords ) / ( RES_ALIGN_SCORE_CONSTANT + rmsd );
				// cerr << "  - Comparing " << pdb_ctr_1;
				// cerr << "\tand " << pdb_ctr_2;
				// cerr << "\t--- num_common_coords : " << num_common_coords;
				// cerr << "\t, rmsd : " << rmsd;
				// cerr << "\t, score : " << score;
				// cerr << "\n";
				edges.emplace_back( pdb_ctr_1, pdb_ctr_2, score );
			}
		}
	}

//	cerr << "Alignment built using residue names is : \n";
//	write_alignment_as_fasta_alignment( cerr, new_alignment, proteins_of_pdbs) << endl;
//	cerr << endl;

	// Return the results
	return make_pair(
		new_alignment,
		get_edges_of_spanning_tree( calc_max_spanning_tree( edges, num_pdbs ) )
	);
}

