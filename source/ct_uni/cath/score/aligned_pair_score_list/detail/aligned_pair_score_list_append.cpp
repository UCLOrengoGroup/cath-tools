/// \file
/// \brief The aligned_pair_score_list_append definitions

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

#include "aligned_pair_score_list_append.hpp"
#include "cath/score/aligned_pair_score/drmsd_score.hpp"
#include "cath/score/aligned_pair_score/gsas_score.hpp"
#include "cath/score/aligned_pair_score/lddt_score.hpp"
#include "cath/score/aligned_pair_score/length_score.hpp"
#include "cath/score/aligned_pair_score/mi_score.hpp"
#include "cath/score/aligned_pair_score/rmsd_score.hpp"
#include "cath/score/aligned_pair_score/sas_score.hpp"
#include "cath/score/aligned_pair_score/sequence_similarity_score.hpp"
#include "cath/score/aligned_pair_score/si_score.hpp"
#include "cath/score/aligned_pair_score/ssap_score.hpp"
#include "cath/score/aligned_pair_score/structal_score.hpp"
#include "cath/score/aligned_pair_score/tm_score.hpp"
#include "cath/structure/protein/amino_acid.hpp"

using namespace cath::score;
using namespace cath::score::detail;

using boost::ptr_vector;

template <>
void score_variety_factory::append_all_varieties<drmsd_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                              ) {
	atom_and_res_pol_varieties_append<drmsd_score>(prm_scores);
}

template <>
void score_variety_factory::append_all_varieties<gsas_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                             ) {
	atom_and_res_pol_varieties_append<gsas_score>(prm_scores);
}

template <>
void score_variety_factory::append_all_varieties<length_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                               ) {
	length_varieties_append<length_score>(prm_scores);
}

template <>
void score_variety_factory::append_all_varieties<lddt_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                             ) {
	threshold_and_atom_and_res_pol_varieties_append<lddt_score>(prm_scores);
}

template <>
void score_variety_factory::append_all_varieties<mi_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                           ) {
	sym_leng_choice_and_atom_and_res_pol_varieties_append<mi_score>(prm_scores);
}

template <>
void score_variety_factory::append_all_varieties<rmsd_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                             ) {
	atom_and_res_pol_varieties_append<rmsd_score>(prm_scores);
}

template <>
void score_variety_factory::append_all_varieties<sequence_similarity_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                                            ) {
	leng_choice_and_subs_matrices_append<sequence_similarity_score>(prm_scores);
}

template <>
void score_variety_factory::append_all_varieties<sas_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                            ) {
	atom_and_res_pol_varieties_append<sas_score>(prm_scores);
}

template <>
void score_variety_factory::append_all_varieties<si_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                           ) {
	sym_leng_choice_and_atom_and_res_pol_varieties_append<si_score>(prm_scores);
}

/// \brief This fails to test many of the varieties of SSAP score so it should be extended
template <>
void score_variety_factory::append_all_varieties<ssap_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                             ) {
	ssap_varieties_append<ssap_score>(prm_scores);
}

/// \brief TODOCUMENT
template <>
void score_variety_factory::append_all_varieties<structal_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                             ) {
	atom_and_res_pol_varieties_append<structal_score>( prm_scores );
}

/// \brief TODOCUMENT
template <>
void score_variety_factory::append_all_varieties<tm_score>(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                           ) {
	atom_and_res_pol_varieties_append<tm_score>( prm_scores );
}

/// Comparing areas under ROC curves on a sensible data-set, an experiment showed that most of the substitution matrices
/// work best when normalised by the number of aligned residues. There were four exceptions:
///  * pam460   got an AUC of 0.709172 with shorter_protein_length        but 0.709163 with num_aligned_residues
///  * pam490   got an AUC of 0.70954  with geometric_mean_protein_length but 0.708105 with num_aligned_residues
///  * pam500   got an AUC of 0.713498 with geometric_mean_protein_length but 0.710388 with num_aligned_residues
///  * identity got an AUC of 0.769681 with longer_protein_length         but 0.704092 with num_aligned_residues
///
/// Apart from "identity", the improvements to be had over num_aligned_residues are very small.
/// Hence this adds the num_aligned_residues version of everything and then adds an extra, standard
/// identity-over-longer_protein_length for good measure.
void score_variety_factory::append_sensible_sequence_similarity_varieties(ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
                                                                          ) {
	subs_matrices_append<sequence_similarity_score>(prm_scores, num_aligned_length_getter() );
	boost::assign::ptr_push_back<sequence_similarity_score>(prm_scores)();
}
