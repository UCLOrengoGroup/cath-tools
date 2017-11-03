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

#include "align_based_superposition_acquirer.hpp"

#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/range/numeric.hpp>
#include <boost/tuple/tuple.hpp>

#include "acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "alignment/alignment.hpp"
#include "alignment/alignment_coord_extractor.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.hpp"
#include "alignment/io/alignment_io.hpp"
#include "alignment/pair_alignment.hpp"
#include "chopping/region/region.hpp"
#include "common/algorithm/sort_copy.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/file/open_fstream.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_list.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "structure/geometry/coord_list.hpp"
#include "superposition/superpose_orient.hpp"
#include "superposition/superposition.hpp"
#include "superposition/superposition_context.hpp"

#include <fstream>
#include <tuple>

using namespace cath;
using namespace cath::align;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::opts;
using namespace cath::sup;

using boost::accumulate;
using boost::filesystem::path;
using std::endl;
using std::get;
using std::ifstream;
using std::make_pair;
using std::make_tuple;
using std::map;
using std::ostream;
using std::pair;
using std::string;
using std::tuple;
using std::vector;

/// \brief TODOCUMENT
superposition_context align_based_superposition_acquirer::do_get_superposition(ostream &arg_stderr ///< TODOCMENT
                                                                               ) const {
	const pdb_list restricted_pdbs = get_restricted_pdbs( *this );
	// Loop over all the PDBs after the first one and grab the common coords between it and the one before
	// and push these onto the data structure that's required for making a superposition
	vector <superposition::indices_and_coord_lists_type> indices_and_coord_lists;
	indices_and_coord_lists.reserve( get_spanning_tree().size() );
	for (const size_size_pair &tree_edge : get_spanning_tree() ) {
		const size_t &index_1 = tree_edge.first;
		const size_t &index_2 = tree_edge.second;

//		arg_stderr << "Extracting common coords between " << name_1 << " and " << name_2 << endl;
//		arg_stderr << "the_alignment.get().num_entries() is " << the_alignment_ref.get().num_entries() << endl;
//		arg_stderr << "the_alignment.get().length()      is " << the_alignment_ref.get().length()      << endl;
//		arg_stderr << "restricted_pdbs.size()            is " << restricted_pdbs.size()                << endl;
//		arg_stderr << "index_1                           is " << index_1                               << endl;
//		arg_stderr << "index_2                           is " << index_2                               << endl;
		const pair<coord_list, coord_list> all_common_coords = alignment_coord_extractor::get_common_coords(
			get_alignment(),
			restricted_pdbs[ index_1 ],
			restricted_pdbs[ index_2 ],
			common_residue_select_all_policy(),
			common_atom_select_ca_policy(),
			index_1,
			index_2
		);
//		const double standard_rmsd_of_original_posns = calc_rmsd(      all_common_coords.first, all_common_coords.second);
		const double standard_rmsd = calc_pairwise_superposition_rmsd( all_common_coords.first, all_common_coords.second);
//		size_t num_within_three_angstroms = 0;
//		size_t num_within_five_angstroms  = 0;
//		for (const size_t &coord_ctr : indices( all_common_coords.first.size() ) ) {
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
			get_alignment(),
			restricted_pdbs,
			index_1,
			index_2
		);

		// std::cerr << "common_coords.first  is " << common_coords.first  << "\n";
		// std::cerr << "cog.first            is " << centre_of_gravity( common_coords.first ) << "\n";
		// std::cerr << "common_coords.second is " << common_coords.second << "\n";

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

		// For reference, epsilon_difference() / relative_difference() on x86_64 Linux 4.8.0-59-generic #64-Ubuntu SMP is :
		//  * double : 2^52
		//  * float  : 2^23
		if ( boost::math::relative_difference( actual_full_rmsd, standard_rmsd ) > TOLERANCE_FOR_EQUAL_RMSDS ) {
			arg_stderr << "Superposed using " << the_selection_policy_acquirer.get_descriptive_name() << " and actual full RMSD is : " << actual_full_rmsd << endl;
		}
	}

	// Construct a superposition and use it to output the PDBs
	const superposition the_superposition = orient_superposition_copy(
		superposition{ indices_and_coord_lists },
		get_alignment(),
		restricted_pdbs
	);

	return {
		the_superposition,
		get_pdbs     ( *this ),
		get_name_sets( *this ),
		get_regions  ( *this ),
		get_alignment(       )
	};
}

/// \brief Ctor for align_based_superposition_acquirer
align_based_superposition_acquirer::align_based_superposition_acquirer(const alignment           &arg_alignment,                ///< TODOCUMENT
                                                                       const size_size_pair_vec  &arg_spanning_tree,            ///< TODOCUMENT
                                                                       const strucs_context      &arg_strucs_context,           ///< TODOCUMENT
                                                                       selection_policy_acquirer  arg_selection_policy_acquirer ///< TODOCUMENT
                                                                       ) : the_alignment_ref             { arg_alignment                              },
                                                                           spanning_tree_ref             { arg_spanning_tree                          },
                                                                           context_ref                   { arg_strucs_context                         },
                                                                           the_selection_policy_acquirer { std::move( arg_selection_policy_acquirer ) } {
}

/// \brief TODOCUMENT
const align::alignment & align_based_superposition_acquirer::get_alignment() const {
	return the_alignment_ref.get();
}

/// \brief TODOCUMENT
const size_size_pair_vec & align_based_superposition_acquirer::get_spanning_tree() const {
	return spanning_tree_ref.get();
}

/// \brief TODOCUMENT
const strucs_context & align_based_superposition_acquirer::get_strucs_context() const {
	return context_ref.get();
}

/// \brief TODOCUMENT
const selection_policy_acquirer & align_based_superposition_acquirer::get_selection_policy_acquirer() const {
	return the_selection_policy_acquirer;
}

/// \brief TODOCUMENT
///
/// \relates align_based_superposition_acquirer
const pdb_list & cath::opts::get_pdbs(const align_based_superposition_acquirer &arg_align_based_superposition_acquirer ///< TODOCUMENT
                                      ) {
	return arg_align_based_superposition_acquirer.get_strucs_context().get_pdbs();
}

/// \brief TODOCUMENT
///
/// \relates align_based_superposition_acquirer
const name_set_list & cath::opts::get_name_sets(const align_based_superposition_acquirer &arg_align_based_superposition_acquirer ///< TODOCUMENT
                                                ) {
	return arg_align_based_superposition_acquirer.get_strucs_context().get_name_sets();
}

/// \brief Getter for the specification of the regions of the PDBs to which the alignment refers
///
/// \relates align_based_superposition_acquirer
const chop::region_vec_opt_vec & cath::opts::get_regions(const align_based_superposition_acquirer &arg_align_based_superposition_acquirer ///< TODOCUMENT
                                                         ) {
	return arg_align_based_superposition_acquirer.get_strucs_context().get_regions();
}

/// \brief Get a copy of the PDBs in the specified align_based_superposition_acquirer, restricted to its regions
///
/// \relates align_based_superposition_acquirer
pdb_list cath::opts::get_restricted_pdbs(const align_based_superposition_acquirer &arg_align_based_superposition_acquirer ///< The align_based_superposition_acquirer to query
                                         ) {
	return get_restricted_pdbs( arg_align_based_superposition_acquirer.get_strucs_context() );
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
	// Store the indices within alignments that each entry/index pair is aligned with something else
	//
	// This is pretty messy but will hopefully not be around too much longer
	//
	// \todo Get rid of hacky_multi_ssap_fuction() and hence rid of this
	using size_size_pair_size_vec_map = map<size_size_pair, size_vec>;
	size_size_pair_size_vec_map aln_idcs_of_entry_and_idx;


	vector <superposition::indices_and_coord_lists_type> indices_and_coord_lists;
	indices_and_coord_lists.reserve(arg_spanning_tree.size());
	for (const size_size_pair &tree_edge : arg_spanning_tree) {
		const bool    flip_names = ( arg_names[ tree_edge.first ] > arg_names[ tree_edge.second ] );
		const size_t &index_1    = flip_names ? tree_edge.second : tree_edge.first;
		const size_t &index_2    = flip_names ? tree_edge.first  : tree_edge.second;
		const string &name_1     = arg_names[ index_1 ];
		const string &name_2     = arg_names[ index_2 ];
		if ( name_1.empty() ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("No name available for " + std::to_string( index_1 ) ));
		}
		if ( name_2.empty() ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("No name available for " + std::to_string( index_2 ) ));
		}

		const path ssap_aln_filename(arg_ssap_align_dir / (name_1 + name_2 + ".list"));
		arg_stderr << "Loading SSAP alignment between " << name_1 << " and " << name_2 << " from " << ssap_aln_filename << endl;
		ifstream ssap_aln_stream;
		open_ifstream(ssap_aln_stream, ssap_aln_filename);
		const alignment the_alignment = read_alignment_from_cath_ssap_legacy_format(
			ssap_aln_stream,
			arg_pdbs[index_1],
			arg_pdbs[index_2],
			ostream_ref{ arg_stderr }
		);
		ssap_aln_stream.close();
//		arg_stderr << "New alignment is " << the_alignment << endl;

//		arg_stderr << "Extracting common coords between " << name_1 << " and " << name_2 << endl;
		check_alignment_is_a_pair( the_alignment );
		for (const auto &aln_index : indices( the_alignment.length() ) ) {
			if ( has_both_positions_of_index( the_alignment, aln_index ) ) {
				const auto entry_indx_a = make_pair( index_1, get_a_position_of_index( the_alignment, aln_index ) );
				aln_idcs_of_entry_and_idx[ entry_indx_a ].push_back( aln_index );
				const auto entry_indx_b = make_pair( index_2, get_b_position_of_index( the_alignment, aln_index ) );
				aln_idcs_of_entry_and_idx[ entry_indx_b ].push_back( aln_index );
			}
		}

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
		// For reference, epsilon_difference() / relative_difference() on x86_64 Linux 4.8.0-59-generic #64-Ubuntu SMP is :
		//  * double : 2^52
		//  * float  : 2^23
		if ( boost::math::relative_difference( actual_full_rmsd, standard_rmsd ) > superposition_acquirer::TOLERANCE_FOR_EQUAL_RMSDS ) {
			arg_stderr << "Superposed using " << arg_selection_policy_acquirer.get_descriptive_name() << " and actual full RMSD is : " << actual_full_rmsd << endl;
		}
	}

	if ( aln_idcs_of_entry_and_idx.empty() ) {
		for (const size_t &entry : indices( arg_pdbs.size() ) ) {
			const pdb &the_pdb = arg_pdbs[ entry ];
			for (const size_t &index : indices( the_pdb.get_num_residues() ) ) {
				aln_idcs_of_entry_and_idx[ make_pair( entry, index ) ].push_back( index );
			}
		}
	}

	const auto mean_aln_idx_entry_index_vec = sort_copy(
		transform_build< vector< tuple< double, size_t, size_t > > >(
			aln_idcs_of_entry_and_idx,
			[] (const auto &x) {
				const auto &entry       = x.first.first;
				const auto &index       = x.first.second;
				const auto &aln_indices = x.second;
				if ( aln_indices.empty() ) {
					BOOST_THROW_EXCEPTION(out_of_range_exception("Cannot find the mean of an empty list of alignment indices"));
				}
				return make_tuple(
					accumulate(
						aln_indices,
						0.0,
						[] (const double &result, const size_t &aln_index) {
							return result + static_cast<double>( aln_index );
						}
					) / static_cast<double>( aln_indices.size() ),
					entry,
					index
				);
			}
		)
	);

	// Construct the unoriented superposition
	const superposition unoriented_supn{ indices_and_coord_lists };

	// Construct an oriented superposition
	return orient_superposition_copy(
		unoriented_supn,
		coord_list{ transform_build<coord_vec>(
			mean_aln_idx_entry_index_vec,
			[&] (const tuple<double, size_t, size_t> &x) {
				return transform_copy(
					unoriented_supn,
					get<1>( x ),
					arg_pdbs[ get<1>( x ) ].get_residue_ca_coord_of_index__backbone_unchecked(
						get<2>( x )
					)
				);
			}
		) }
	);
}
