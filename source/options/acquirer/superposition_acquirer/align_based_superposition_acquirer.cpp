/// \file
/// \brief The align_based_superposition_acquirer class definitions

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

#include "align_based_superposition_acquirer.h"

#include <boost/lexical_cast.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/tuple/tuple.hpp>

#include "alignment/alignment.h"
#include "alignment/alignment_coord_extractor.h"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.h"
#include "alignment/io/alignment_io.h"
#include "common/file/open_fstream.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"
#include "options/acquirer/alignment_acquirer/alignment_acquirer.h"
#include "options/executable/cath_superpose_options/cath_superpose_options.h"
#include "structure/geometry/coord_list.h"
#include "superposition/superposition.h"
#include "superposition/superposition_context.h"

#include <fstream>
#include <iostream> // ***** TEMPORARY *****

using namespace boost::filesystem;
using namespace boost::test_tools;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

using boost::lexical_cast;

const pdb_list & align_based_superposition_acquirer::get_pdbs_cref() const {
	return pdbs;
}

const str_vec & align_based_superposition_acquirer::get_names_cref() const {
	return names;
}

/// \brief TODOCUMENT
superposition_context align_based_superposition_acquirer::do_get_superposition(ostream &arg_stderr
                                                                               ) const {
	// Loop over all the PDBs after the first one and grab the common coords between it and the one before
	// and push these onto the data structure that's required for making a superposition
	vector <superposition::indices_and_coord_lists_type> indices_and_coord_lists;
	indices_and_coord_lists.reserve(spanning_tree.size());
	for (const size_size_pair &tree_edge : spanning_tree) {
		const size_t &index_1 = tree_edge.first;
		const size_t &index_2 = tree_edge.second;
		const string &name_1  = names[ index_1 ];
		const string &name_2  = names[ index_2 ];
		if (name_1.empty()) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("No name available for " + lexical_cast<string>( index_1 ) ));
		}
		if (name_2.empty()) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("No name available for " + lexical_cast<string>( index_2 ) ));
		}

//		arg_stderr << "Extracting common coords between " << name_1 << " and " << name_2 << endl;
//		arg_stderr << "the_alignment.num_entries() is " << the_alignment.num_entries()   << endl;
//		arg_stderr << "the_alignment.length()      is " << the_alignment.length()        << endl;
//		arg_stderr << "pdbs.size()                 is " << pdbs.size()                   << endl;
//		arg_stderr << "index_1                     is " << index_1                       << endl;
//		arg_stderr << "index_2                     is " << index_2                       << endl;
		const pair<coord_list, coord_list> all_common_coords = alignment_coord_extractor::get_common_coords(
			the_alignment,
			pdbs[ index_1 ],
			pdbs[ index_2 ],
			common_residue_select_all_policy(),
			common_atom_select_ca_policy(),
			index_1,
			index_2
		);
//		const double standard_rmsd_of_original_posns = calc_rmsd(      all_common_coords.first, all_common_coords.second);
		const double standard_rmsd = calc_pairwise_superposition_rmsd( all_common_coords.first, all_common_coords.second);
//		size_t num_within_three_angstroms = 0;
//		size_t num_within_five_angstroms  = 0;
//		for (size_t coord_ctr = 0; coord_ctr < all_common_coords.first.size(); ++coord_ctr) {
//			const coord coord_a = all_common_coords.first[coord_ctr];
//			const coord coord_b = all_common_coords.second[coord_ctr];
//			const double distance = distance_between_points(coord_a, coord_b);
////			arg_stderr << "distance : " << distance << endl;
//			if (distance <= 3.0) {
//				++num_within_three_angstroms;
//			}
//			if (distance <= 5.0) {
//				++num_within_five_angstroms;
//			}
//		}

//		arg_stderr << "Standard RMSD of the structures in their original positions is : " << standard_rmsd_of_original_posns << endl;
		arg_stderr << "Standard RMSD is : " << standard_rmsd << endl;
//		arg_stderr << "percentage_within_three_angstrom       " << (100.0 * num_within_three_angstroms  / all_common_coords.first.size()) << endl;;
//		arg_stderr << " (" << (100.0 * num_within_three_angstroms / all_common_coords.first.size()) << " %)" << endl;
//		arg_stderr << "percentage_within_five_angstroms " << (100.0 * num_within_five_angstroms  / all_common_coords.first.size()) << endl;;
//		arg_stderr << " (" << (100.0 * num_within_five_angstroms  / all_common_coords.first.size()) << " %)" << endl;

//		const common_residue_selection_policy &policy = (
//			false
////			the_alignment.is_scored()
//			? static_cast<const common_residue_selection_policy &>(best_score_percent_policy)
//			: static_cast<const common_residue_selection_policy &>(select_all_policy)
//		);

//		const selection_policy_acquirer the_selection_policy_acquirer = arg_cath_superpose_options.get_selection_policy_acquirer();
		const pair<coord_list, coord_list> common_coords = get_common_coords(
			the_selection_policy_acquirer,
			the_alignment,
			pdbs,
			index_1,
			index_2
		);

		// Add the data to the back of the data to be used to make a superposition
		indices_and_coord_lists.emplace_back(
			index_1,
			common_coords.first,
			index_2,
			common_coords.second
		);

		const superposition pairwise_sup = create_pairwise_superposition(common_coords.first, common_coords.second);
		const double actual_full_rmsd = calc_rmsd_between_superposed_entries(
			pairwise_sup,
			superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION,  all_common_coords.first,
			superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION, all_common_coords.second
		);
		if (!check_is_close(actual_full_rmsd, standard_rmsd, percent_tolerance(PERCENT_TOLERANCE_FOR_EQUAL_RMSDS))) {
			arg_stderr << "Superposed using " << the_selection_policy_acquirer.get_descriptive_name() << " and actual full RMSD is : " << actual_full_rmsd << endl;
		}
	}

	// Construct a superposition and use it to output the PDBs
	const superposition the_superposition(indices_and_coord_lists);

	return superposition_context(pdbs, names, the_superposition, the_alignment);
}

/// \brief Ctor for align_based_superposition_acquirer
align_based_superposition_acquirer::align_based_superposition_acquirer(const pdb_list                  &arg_pdbs,                     ///< TODOCUMENT
                                                                       const str_vec                   &arg_names,                    ///< TODOCUMENT
                                                                       const alignment                 &arg_alignment,                ///< TODOCUMENT
                                                                       const size_size_pair_vec        &arg_spanning_tree,            ///< TODOCUMENT
                                                                       const selection_policy_acquirer &arg_selection_policy_acquirer ///< TODOCUMENT
                                                                       ) : pdbs(arg_pdbs),
                                                                           names(arg_names),
                                                                           the_alignment(arg_alignment),
                                                                           spanning_tree(arg_spanning_tree),
                                                                           the_selection_policy_acquirer(arg_selection_policy_acquirer) {
}

/// \relates align_based_superposition_acquirer
///
/// REALLY HACKY CODE TO EXPLORE MULTIPLE SSAP SUPERPOSITIONS
/// \todo REMOVE THIS AND IMPLEMENT IT PROPERLY IF IT'S ANY GOOD
///       STEPS:
///        - write something to glue alignments together
///        - might need something to exchange positions of alignments?
superposition cath::opts::hacky_multi_ssap_fuction(const pdb_list                  &arg_pdbs,                      ///< TODOCUMENT
                                                   const str_vec                   &arg_names,                     ///< TODOCUMENT
                                                   const size_size_pair_vec        &arg_spanning_tree,             ///< TODOCUMENT
                                                   const path                      &arg_ssap_align_dir,            ///< TODOCUMENT
                                                   const selection_policy_acquirer &arg_selection_policy_acquirer, ///< TODOCUMENT
                                                   ostream                         &arg_stderr                     ///< TODOCUMENT
                                                   ) {
	vector <superposition::indices_and_coord_lists_type> indices_and_coord_lists;
	indices_and_coord_lists.reserve(arg_spanning_tree.size());
	for (const size_size_pair &tree_edge : arg_spanning_tree) {
		const bool    flip_names = ( arg_names[ tree_edge.first ] > arg_names[ tree_edge.second ] );
		const size_t &index_1    = flip_names ? tree_edge.second : tree_edge.first;
		const size_t &index_2    = flip_names ? tree_edge.first  : tree_edge.second;
		const string &name_1     = arg_names[ index_1 ];
		const string &name_2     = arg_names[ index_2 ];
		if ( name_1.empty() ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("No name available for " + lexical_cast<string>( index_1 ) ));
		}
		if ( name_2.empty() ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("No name available for " + lexical_cast<string>( index_2 ) ));
		}

		const path ssap_aln_filename(arg_ssap_align_dir / (name_1 + name_2 + ".list"));
		arg_stderr << "Loading SSAP alignment between " << name_1 << " and " << name_2 << " from " << ssap_aln_filename << endl;
		ifstream ssap_aln_stream;
		open_ifstream(ssap_aln_stream, ssap_aln_filename);
		const alignment the_alignment = read_alignment_from_cath_ssap_legacy_format(
			ssap_aln_stream,
			arg_pdbs[index_1],
			arg_pdbs[index_2],
			arg_stderr
		);
		ssap_aln_stream.close();
//		arg_stderr << "New alignment is " << the_alignment << endl;

//		arg_stderr << "Extracting common coords between " << name_1 << " and " << name_2 << endl;
		const common_residue_select_best_score_percent_policy best_score_percent_policy;
		const common_residue_select_all_policy                select_all_policy{};
		const pair<coord_list, coord_list> all_common_coords = alignment_coord_extractor::get_common_coords(
			the_alignment,
			arg_pdbs[ index_1 ],
			arg_pdbs[ index_2 ],
			select_all_policy,
			common_atom_select_ca_policy(),
			alignment::PAIR_A_IDX,
			alignment::PAIR_B_IDX
		);
		const double standard_rmsd = calc_pairwise_superposition_rmsd(all_common_coords.first, all_common_coords.second);
		arg_stderr << "Standard RMSD is : " << standard_rmsd << endl;

		const pair<coord_list, coord_list> common_coords = arg_selection_policy_acquirer.get_common_coords(
			the_alignment,
			arg_pdbs[ index_1 ],
			arg_pdbs[ index_2 ],
			alignment::PAIR_A_IDX,
			alignment::PAIR_B_IDX
		);

		// Add the data to the back of the data to be used to make a superposition
		indices_and_coord_lists.emplace_back(
			index_1,
			common_coords.first,
			index_2,
			common_coords.second
		);

		const superposition pairwise_sup = create_pairwise_superposition(common_coords.first, common_coords.second);
		const double actual_full_rmsd = calc_rmsd_between_superposed_entries(
			pairwise_sup,
			superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION,  all_common_coords.first,
			superposition::INDEX_OF_SECOND_IN_PAIRWISE_SUPERPOSITION, all_common_coords.second
		);
		if (!check_is_close(actual_full_rmsd, standard_rmsd, percent_tolerance(superposition_acquirer::PERCENT_TOLERANCE_FOR_EQUAL_RMSDS))) {
			arg_stderr << "Superposed using " << arg_selection_policy_acquirer.get_descriptive_name() << " and actual full RMSD is : " << actual_full_rmsd << endl;
		}
	}
	return superposition(indices_and_coord_lists);
}
