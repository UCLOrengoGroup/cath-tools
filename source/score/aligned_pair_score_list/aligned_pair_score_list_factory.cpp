/// \file
/// \brief The aligned_pair_score_list_factory class definitions

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

#define BOOST_ASSIGN_MAX_PARAMS 7

#include "aligned_pair_score_list_factory.hpp"

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/common_atom_selection_policy/common_atom_select_ca_policy.hpp"
#include "alignment/common_atom_selection_policy/common_atom_select_cb_policy.hpp"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_min_score_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "exception/not_implemented_exception.hpp"
#include "score/aligned_pair_score/sequence_similarity_score.hpp"
#include "score/aligned_pair_score/drmsd_score.hpp"
#include "score/aligned_pair_score/gsas_score.hpp"
#include "score/aligned_pair_score/mi_score.hpp"
#include "score/aligned_pair_score/lddt_score.hpp"
#include "score/aligned_pair_score/length_score.hpp"
#include "score/aligned_pair_score/overlap_score.hpp"
#include "score/aligned_pair_score/rmsd_score.hpp"
#include "score/aligned_pair_score/sas_score.hpp"
#include "score/aligned_pair_score/si_score.hpp"
#include "score/aligned_pair_score/ssap_score.hpp"
#include "score/aligned_pair_score/structal_score.hpp"
#include "score/aligned_pair_score/substitution_matrix/blosum62_substitution_matrix.hpp"
#include "score/aligned_pair_score/tm_score.hpp"
#include "score/aligned_pair_score_list/aligned_pair_score_list.hpp"
#include "score/aligned_pair_score_list/detail/aligned_pair_score_list_append.hpp"
#include "score/length_getter/geometric_mean_length_getter.hpp"
#include "score/length_getter/length_of_first_getter.hpp"
#include "score/length_getter/length_of_longer_getter.hpp"
#include "score/length_getter/length_of_second_getter.hpp"
#include "score/length_getter/length_of_shorter_getter.hpp"
#include "score/length_getter/mean_length_getter.hpp"
#include "structure/protein/amino_acid.hpp"

#include <vector>

using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace cath::score::detail;
using namespace std;

using boost::assign::ptr_push_back;
using boost::ptr_vector;

/// \brief Factory function for making an aligned_pair_score_list that contains every type of aligned_pair_score there is
///
/// \todo This originally tried to use boost::mpl::for_each for each group of aligned_pair_scores
///       that take the same possible parameters but it got a bit ugly and wasn't worth the effort.
///       Still, there's almost certainly a nicer way to do this, along those lines.
///
/// \relates aligned_pair_score_list
aligned_pair_score_list cath::score::make_full_aligned_pair_score_list() {
	// Create a ptr_vector of scores to be populated
	ptr_vector<aligned_pair_score> scores;

	score_variety_factory::append_all_varieties< drmsd_score               >( scores );
	score_variety_factory::append_all_varieties< gsas_score                >( scores );
	score_variety_factory::append_all_varieties< length_score              >( scores );
	score_variety_factory::append_all_varieties< lddt_score                >( scores );
	score_variety_factory::append_all_varieties< mi_score                  >( scores );
	score_variety_factory::append_all_varieties< ssap_score                >( scores );
	score_variety_factory::append_all_varieties< rmsd_score                >( scores );
	score_variety_factory::append_all_varieties< sas_score                 >( scores );
	score_variety_factory::append_all_varieties< si_score                  >( scores );

	// ***** TEMPORARILY REPLACE THE FULL LIST OF SEQUENCE_SIMILARITY_SCORES BLOSUM62 AND IDENTITY
//	score_variety_factory::append_all_varieties< sequence_similarity_score >( scores );
	ptr_push_back< sequence_similarity_score >( scores )( length_of_longer_getter(),  make_subs_matrix_blosum62() );
	ptr_push_back< sequence_similarity_score >( scores )( length_of_longer_getter(),  make_subs_matrix_identity() );

	score_variety_factory::append_all_varieties< structal_score            >( scores );
	score_variety_factory::append_all_varieties< tm_score                  >( scores );

	ptr_push_back< naive_overlap  >( scores )();
	ptr_push_back< local_overlap  >( scores )();
	ptr_push_back< global_overlap >( scores )();

	// Return an aligned_pair_score_list of the scores
	const aligned_pair_score_list results( scores );
	warn_on_duplicate_human_friendly_names( results );
	return results;
}

aligned_pair_score_list cath::score::make_seq_sim_flavours_aligned_pair_score_list() {
	// Create a ptr_vector of scores to be populated
	ptr_vector<aligned_pair_score> scores;

	score_variety_factory::append_all_varieties< sequence_similarity_score >( scores );

	// Return an aligned_pair_score_list of the scores
	const aligned_pair_score_list results( scores );
	warn_on_duplicate_human_friendly_names( results );
	return results;
}

aligned_pair_score_list cath::score::make_ssap_flavours_aligned_pair_score_list() {
	// Create a ptr_vector of scores to be populated
	ptr_vector<aligned_pair_score> scores;

	score_variety_factory::append_sensible_sequence_similarity_varieties( scores );

	// Return an aligned_pair_score_list of the scores
	const aligned_pair_score_list results( scores );
	warn_on_duplicate_human_friendly_names( results );
	return results;
}

/// \brief Factory function for making an aligned_pair_score_list that contains one entry of each concrete type
///
/// \relates aligned_pair_score_list
aligned_pair_score_list cath::score::make_one_of_each_type_aligned_pair_score_list() {
	// Create a ptr_vector of scores to be populated
	ptr_vector<aligned_pair_score> scores;

	score_variety_factory::append_all_varieties< length_score >( scores );

	ptr_push_back< drmsd_score               >( scores )(                           );
	ptr_push_back< gsas_score                >( scores )(                           );
	ptr_push_back< lddt_score                >( scores )(                           );
	ptr_push_back< mi_score                  >( scores )(                           ); // MI
	ptr_push_back< mi_score                  >( scores )( length_of_longer_getter() ); // MIMAX
	ptr_push_back< ssap_score                >( scores )(                           );
	ptr_push_back< rmsd_score                >( scores )(                           );
	ptr_push_back< sas_score                 >( scores )(                           );
	ptr_push_back< si_score                  >( scores )(                           ); // SI
	ptr_push_back< si_score                  >( scores )( length_of_longer_getter() ); // SIMAX
	ptr_push_back< sequence_similarity_score >( scores )(                           );
	ptr_push_back< tm_score                  >( scores )(                           );
	ptr_push_back< naive_overlap             >( scores )(                           );
	ptr_push_back< local_overlap             >( scores )(                           );
	ptr_push_back< global_overlap            >( scores )(                           );

	score_variety_factory::append_sensible_sequence_similarity_varieties( scores );

	// Return an aligned_pair_score_list of the scores
	const aligned_pair_score_list results( scores );
	warn_on_duplicate_human_friendly_names( results );
	return results;
}


/// \brief Factory function for making an aligned_pair_score_list that contains a sensible list of defaults
///
/// \relates aligned_pair_score_list
aligned_pair_score_list cath::score::make_default_aligned_pair_score_list() {
	const common_residue_select_all_policy select_all_residues{};
//	const common_atom_select_ca_policy     select_ca_atoms;
	const common_atom_select_cb_policy     select_cb_atoms{};

	// Create a ptr_vector of scores to be populated
	ptr_vector<aligned_pair_score> scores;

	ptr_push_back< length_score              >( scores )( length_of_longer_getter()                                                                                    ); // AUC:0.598147 longer_protein_length
	ptr_push_back< length_score              >( scores )( length_of_shorter_getter()                                                                                   ); // AUC:0.696794 shorter_protein_length
	ptr_push_back< length_score              >( scores )( num_aligned_length_getter()                                                                                  ); // AUC:0.721526 num_aligned_residues

	ptr_push_back< length_score              >( scores )( geometric_mean_length_getter()                                                                               ); // AUC:0.663327 geometric_mean_protein_length
	ptr_push_back< length_score              >( scores )( mean_length_getter()                                                                                         ); // AUC:0.646203 mean_protein_length

	ptr_push_back< ssap_score                >( scores )(                                                                                                              ); // AUC:0.815189 ssap
	ptr_push_back< drmsd_score               >( scores )(                                                                                                              ); // AUC:0.706255 dRMSD
	ptr_push_back< rmsd_score                >( scores )(                                                                                                              ); // AUC:0.703609 RMSD
	ptr_push_back< gsas_score                >( scores )(                                                                                                              ); // AUC:0.798789 SAS
	ptr_push_back< sas_score                 >( scores )(                                                                                                              ); // AUC:0.800445 GSAS
	ptr_push_back< si_score                  >( scores )( length_of_longer_getter()                                                                                    ); // AUC:0.767007 SIMAX
	ptr_push_back< si_score                  >( scores )( length_of_shorter_getter()                                                                                   ); // AUC:0.708524 SI
	ptr_push_back< mi_score                  >( scores )( length_of_longer_getter()                                                                                    ); // AUC:0.769700 MIMAX
	ptr_push_back< mi_score                  >( scores )( length_of_shorter_getter()                                                                                   ); // AUC:0.707264 MI
	ptr_push_back< lddt_score                >( scores )( lddt_distance_threshold::DEFAULT_AVG                                                                         ); // AUC:0.764472 lDDT.threshold[STD_MEAN]
	ptr_push_back< lddt_score                >( scores )( lddt_distance_threshold::HALF_A                                                                              ); // AUC:0.713065 lDDT.threshold[0.5]
	ptr_push_back< lddt_score                >( scores )( lddt_distance_threshold::ONE_A                                                                               ); // AUC:0.750847 lDDT.threshold[1]
	ptr_push_back< lddt_score                >( scores )( lddt_distance_threshold::TWO_A                                                                               ); // AUC:0.773225 lDDT.threshold[2]
	ptr_push_back< lddt_score                >( scores )( lddt_distance_threshold::FOUR_A                                                                              ); // AUC:0.773487 lDDT.threshold[4]
	ptr_push_back< structal_score            >( scores )(                                                                                                              ); // AUC:0.800231 structal
	ptr_push_back< tm_score                  >( scores )(                                                                                                              ); // AUC:0.797982 TM-score
	ptr_push_back< naive_overlap             >( scores )(                                                                                                              ); // AUC:0.640879 overlap.shorter_protein_length.longer_protein_length
	ptr_push_back< local_overlap             >( scores )(                                                                                                              ); // AUC:0.563231 overlap.num_aligned_residues.shorter_protein_length
	ptr_push_back< global_overlap            >( scores )(                                                                                                              ); // AUC:0.699839 overlap.num_aligned_residues.longer_protein_length

	ptr_push_back< sequence_similarity_score >( scores )( length_of_longer_getter(),  make_subs_matrix_blosum62()                                                      ); // AUC:0.633739 sequence_id.blosum62
	ptr_push_back< sequence_similarity_score >( scores )( length_of_longer_getter(),  make_subs_matrix_identity()                                                      ); // AUC:0.771108 sequence_id.identity

	ptr_push_back< si_score                  >( scores )( length_of_longer_getter(), common_residue_select_best_score_percent_policy(), common_atom_select_cb_policy() ); // AUC:0.796550 SIMAX.select_best_score_percent[70].cb_atoms
	ptr_push_back< mi_score                  >( scores )( length_of_longer_getter(), common_residue_select_best_score_percent_policy(), common_atom_select_cb_policy() ); // AUC:0.798375 MIMAX.select_best_score_percent[70].cb_atoms
	ptr_push_back< tm_score                  >( scores )(                            common_residue_select_min_score_policy( 0.01 ),    common_atom_select_ca_policy() ); // AUC:0.798518 TM-score.select_min_score[0.01]
	ptr_push_back< sas_score                 >( scores )(                            common_residue_select_min_score_policy( 0.01 ),    common_atom_select_ca_policy() ); // AUC:0.799165 SAS.select_min_score[0.01]
	ptr_push_back< sas_score                 >( scores )(                            common_residue_select_all_policy(),                common_atom_select_cb_policy() ); // AUC:0.800209 SAS.cb_atoms
	ptr_push_back< sas_score                 >( scores )(                            common_residue_select_min_score_policy( 0.01 ),    common_atom_select_cb_policy() ); // AUC:0.800528 SAS.select_min_score[0.01].cb_atoms
	ptr_push_back< gsas_score                >( scores )(                            common_residue_select_min_score_policy( 0.01 ),    common_atom_select_ca_policy() ); // AUC:0.800802 GSAS.select_min_score[0.01]
	ptr_push_back< structal_score            >( scores )(                            common_residue_select_min_score_policy( 0.01 ),    common_atom_select_ca_policy() ); // AUC:0.801087 structal.select_min_score[0.01]
	ptr_push_back< tm_score                  >( scores )(                            common_residue_select_all_policy(),                common_atom_select_cb_policy() ); // AUC:0.801571 TM-score.cb_atoms
	ptr_push_back< gsas_score                >( scores )(                            common_residue_select_all_policy(),                common_atom_select_cb_policy() ); // AUC:0.801881 GSAS.cb_atoms
	ptr_push_back< tm_score                  >( scores )(                            common_residue_select_min_score_policy( 0.01 ),    common_atom_select_cb_policy() ); // AUC:0.802079 TM-score.select_min_score[0.01].cb_atoms
	ptr_push_back< gsas_score                >( scores )(                            common_residue_select_min_score_policy( 0.01 ),    common_atom_select_cb_policy() ); // AUC:0.802196 GSAS.select_min_score[0.01].cb_atoms
	ptr_push_back< structal_score            >( scores )(                            common_residue_select_all_policy(),                common_atom_select_cb_policy() ); // AUC:0.803037 structal.cb_atoms
	ptr_push_back< structal_score            >( scores )(                            common_residue_select_min_score_policy( 0.01 ),    common_atom_select_cb_policy() ); // AUC:0.803875 structal.select_min_score[0.01].cb_atoms
	ptr_push_back< structal_score            >( scores )(                            common_residue_select_best_score_percent_policy(), common_atom_select_ca_policy() ); // AUC:0.817068 structal.select_best_score_percent[70]
	ptr_push_back< structal_score            >( scores )(                            common_residue_select_best_score_percent_policy(), common_atom_select_cb_policy() ); // AUC:0.818635 structal.select_best_score_percent[70].cb_atoms
	ptr_push_back< sas_score                 >( scores )(                            common_residue_select_best_score_percent_policy(), common_atom_select_ca_policy() ); // AUC:0.822474 SAS.select_best_score_percent[70]
	ptr_push_back< gsas_score                >( scores )(                            common_residue_select_best_score_percent_policy(), common_atom_select_ca_policy() ); // AUC:0.822660 GSAS.select_best_score_percent[70]
	ptr_push_back< sas_score                 >( scores )(                            common_residue_select_best_score_percent_policy(), common_atom_select_cb_policy() ); // AUC:0.824783 SAS.select_best_score_percent[70].cb_atoms
	ptr_push_back< gsas_score                >( scores )(                            common_residue_select_best_score_percent_policy(), common_atom_select_cb_policy() ); // AUC:0.825142 GSAS.select_best_score_percent[70].cb_atoms

	ptr_push_back< ssap_score                >( scores )( length_of_longer_getter(), common_residue_select_all_policy(),                common_atom_select_cb_policy(),
	                                                      ssap_score::default_post_processing, ssap_score_accuracy::HIGH, 10, distance_score_formula::SIMPLIFIED       ); // AUC:0.824031 ssap.cb_atoms.high_accuracy.num_excluded_on_sides:10.distance_score_formula_simplified

	// Return an aligned_pair_score_list of the scores
	const aligned_pair_score_list results( scores );
	warn_on_duplicate_human_friendly_names( results );
	return results;
}

/// \brief Factory function for making an aligned_pair_score_list that contains those aligned_pair_scores used in the old SSAP output
///
/// Note that this does not include the two names that typically appear at the start of the old SSAP scores output
///
/// The scores are:
///  * length1    : Length of protein 1
///  * length2    : Length of protein 2
///  * ssap_score : SSAP score for structural comparison (0-100)
///  * num_equivs : Number of equivalent/aligned residues
///  * overlap_pc : Percentage overlap  (100% x overlap /length of largest)
///  * seq_id_pc  : Percentage identity (100% x identity/length of smallest)
///  * rmsd       : RMSD of superposed structures
///
/// \relates aligned_pair_score_list
aligned_pair_score_list cath::score::make_old_ssap_aligned_pair_score_list() {
	// Create a ptr_vector of scores to be populated
	ptr_vector<aligned_pair_score> scores;

	ptr_push_back< length_score              >( scores )( length_of_first_getter()    ); //  length1    : Length of protein 1
	ptr_push_back< length_score              >( scores )( length_of_second_getter()   ); //  length2    : Length of protein 2
	ptr_push_back< ssap_score                >( scores )(                             ); //  ssap_score : SSAP score for structural comparison (0-100)
	ptr_push_back< length_score              >( scores )( num_aligned_length_getter() ); //  num_equivs : Number of equivalent/aligned residues
	ptr_push_back< global_overlap            >( scores )(                             ); //  overlap_pc : Percentage overlap  (100% x overlap /length of largest)
	ptr_push_back< sequence_similarity_score >( scores )(                             ); //  seq_id_pc  : Percentage identity (100% x identity/length of smallest)
	ptr_push_back< rmsd_score                >( scores )(                             ); //  rmsd       : RMSD of superposed structures

	// Return an aligned_pair_score_list of the scores
	const aligned_pair_score_list results( scores );
	warn_on_duplicate_human_friendly_names( results );
	return results;
}
