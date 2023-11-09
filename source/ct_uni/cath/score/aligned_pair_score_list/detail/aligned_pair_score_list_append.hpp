/// \file
/// \brief The aligned_pair_score_list_append header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_DETAIL_ALIGNED_PAIR_SCORE_LIST_APPEND_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_DETAIL_ALIGNED_PAIR_SCORE_LIST_APPEND_HPP

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "cath/alignment/common_atom_selection_policy/common_atom_select_cb_policy.hpp" // ***** TEMPORARY *****
#include "cath/alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp" // ***** TEMPORARY *****
#include "cath/alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/score/aligned_pair_score/lddt_score.hpp"
#include "cath/score/aligned_pair_score/ssap_score.hpp"
#include "cath/score/aligned_pair_score/ssap_score/ssap_score_post_processing.hpp"
#include "cath/score/aligned_pair_score/substitution_matrix/substitution_matrix.hpp"
#include "cath/score/length_getter/sym_protein_only_length_getter.hpp"
#include "cath/ssap/distance_score_formula.hpp"

#include <iostream> // ***** TEMPORARY *****

// clang-format off
namespace cath::align { class common_atom_selection_policy; }
namespace cath::align { class common_residue_selection_policy; }
namespace cath::score { class aligned_pair_score; }
namespace cath::score { class aligned_pair_score_list; }
namespace cath::score { class drmsd_score; }
namespace cath::score { class gsas_score; }
namespace cath::score { class lddt_score; }
namespace cath::score { class length_getter; }
namespace cath::score { class length_score; }
namespace cath::score { class mi_score; }
namespace cath::score { class rmsd_score; }
namespace cath::score { class sas_score; }
namespace cath::score { class sequence_similarity_score; }
namespace cath::score { class si_score; }
namespace cath::score { class ssap_score; }
namespace cath::score { class structal_score; }
namespace cath::score { class tm_score; }
// clang-format on

namespace cath::score::detail {

	using cath::common::literals::operator""_z;

	/// \brief TODOCUMENT
	class score_variety_factory final {
	private:
		template <typename T>
		static void length_varieties_append(boost::ptr_vector<aligned_pair_score> &);

		template <typename T>
		static void atom_and_res_pol_varieties_append(boost::ptr_vector<aligned_pair_score> &);

		template <typename T>
		static void sym_leng_choice_and_atom_and_res_pol_varieties_append(boost::ptr_vector<aligned_pair_score> &);

		template <typename T>
		static void threshold_and_atom_and_res_pol_varieties_append(boost::ptr_vector<aligned_pair_score> &);

		template <typename T>
		static void ssap_varieties_append(boost::ptr_vector<aligned_pair_score> &);

		template <typename T>
		static void leng_choice_and_subs_matrices_append(boost::ptr_vector<aligned_pair_score> &);

		template <typename T>
		static void subs_matrices_append(boost::ptr_vector<aligned_pair_score> &,
		                                 const length_getter &);

	public:
		score_variety_factory() = delete;

		template <typename T>
		static void append_all_varieties(boost::ptr_vector<aligned_pair_score> &);

		static void append_sensible_sequence_similarity_varieties(boost::ptr_vector<aligned_pair_score> &);
	};

	/// \brief TODOCUMENT
	template <typename T>
	void score_variety_factory::length_varieties_append(boost::ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
	                                                    ) {
		const boost::ptr_vector<length_getter> all_length_getters = get_all_length_getters();
		for (const length_getter &length_getter : all_length_getters) {
			boost::assign::ptr_push_back< T >( prm_scores )( length_getter );
		}
	}

	/// \brief Append a score object of the specified type to the specified ptr_vector of scores
	template <typename T>
	void score_variety_factory::atom_and_res_pol_varieties_append(boost::ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
	                                                              ) {
		const boost::ptr_vector<align::common_residue_selection_policy> all_comm_res_seln_pols  = cath::align::get_all_common_residue_selection_policies();
		const boost::ptr_vector<align::common_atom_selection_policy>    all_comm_atom_seln_pols = cath::align::get_all_common_atom_selection_policies();
		for (const align::common_residue_selection_policy &comm_res_seln_pol : all_comm_res_seln_pols) {
			for (const align::common_atom_selection_policy &comm_atom_seln_pol : all_comm_atom_seln_pols) {
				boost::assign::ptr_push_back< T >( prm_scores )( comm_res_seln_pol, comm_atom_seln_pol );
			}
		}
	}

	/// \brief Append all varieties of score object of the specified type to the specified ptr_vector of scores
	template <typename T>
	void score_variety_factory::sym_leng_choice_and_atom_and_res_pol_varieties_append(boost::ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
	                                                                                  ) {
		const boost::ptr_vector<sym_protein_only_length_getter>         all_sym_length_choices  = cath::score::get_all_sym_protein_only_length_getters();
		const boost::ptr_vector<align::common_residue_selection_policy> all_comm_res_seln_pols  = cath::align::get_all_common_residue_selection_policies();
		const boost::ptr_vector<align::common_atom_selection_policy>    all_comm_atom_seln_pols = cath::align::get_all_common_atom_selection_policies();
		for (const sym_protein_only_length_getter &sym_length_choice : all_sym_length_choices) {
			for (const align::common_residue_selection_policy &comm_res_seln_pol : all_comm_res_seln_pols) {
				for (const align::common_atom_selection_policy &comm_atom_seln_pol : all_comm_atom_seln_pols) {
					boost::assign::ptr_push_back< T >( prm_scores )( sym_length_choice, comm_res_seln_pol, comm_atom_seln_pol );
				}
			}
		}
	}

	/// \brief Append all varieties of score object of the specified type to the specified ptr_vector of scores
	template <typename T>
	void score_variety_factory::threshold_and_atom_and_res_pol_varieties_append(boost::ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
	                                                                            ) {
		const std::vector<lddt_distance_threshold>                      all_lddt_distance_thresholds = cath::score::get_all_lddt_distance_thresholds();
		// std::cerr << "all_lddt_distance_thresholds has size : " << all_lddt_distance_thresholds.size() << std::endl;
		const boost::ptr_vector<align::common_residue_selection_policy> all_comm_res_seln_pols       = cath::align::get_all_common_residue_selection_policies();
		// std::cerr << "all_comm_res_seln_pols       has size : " << all_comm_res_seln_pols.size()       << std::endl;
		const boost::ptr_vector<align::common_atom_selection_policy>    all_comm_atom_seln_pols      = cath::align::get_all_common_atom_selection_policies();
		// std::cerr << "all_comm_atom_seln_pols      has size : " << all_comm_atom_seln_pols.size()      << std::endl;
		for (const lddt_distance_threshold &lddt_dist_thresh : all_lddt_distance_thresholds) {
			for (const align::common_residue_selection_policy &comm_res_seln_pol : all_comm_res_seln_pols) {
				for (const align::common_atom_selection_policy &comm_atom_seln_pol : all_comm_atom_seln_pols) {
					boost::assign::ptr_push_back< T >( prm_scores )( lddt_dist_thresh, comm_res_seln_pol, comm_atom_seln_pol );
				}
			}
		}
	}

	/// \brief Append all varieties of score object of the specified type to the specified ptr_vector of scores
	template <typename T>
	void score_variety_factory::ssap_varieties_append(boost::ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
	                                                  ) {
		const boost::ptr_vector<sym_protein_only_length_getter>         all_sym_length_choices  = cath::score::get_all_sym_protein_only_length_getters();
//		const boost::ptr_vector<align::common_residue_selection_policy> all_comm_res_seln_pols  = cath::align::get_all_common_residue_selection_policies();
//		const boost::ptr_vector<align::common_atom_selection_policy>    all_comm_atom_seln_pols = cath::align::get_all_common_atom_selection_policies();

		boost::assign::ptr_push_back< T >( prm_scores )( );

		for (const sym_protein_only_length_getter &sym_length_choice : all_sym_length_choices) {
//			for (const align::common_residue_selection_policy &comm_res_seln_pol : all_comm_res_seln_pols) {
//				for (const align::common_atom_selection_policy &comm_atom_seln_pol : all_comm_atom_seln_pols) {
					for (const ssap_score_post_processing &post_proc : all_ssap_score_post_processings) {
//						for (const ssap_score_accuracy &accuracy : all_ssap_score_accuracies) {
							for (const size_t num : { 0_z, 5_z, 10_z, 15_z } ) {
								for (const distance_score_formula &dist_score : all_distance_score_formulae) {
									boost::assign::ptr_push_back< T >( prm_scores )(
										sym_length_choice,
										align::common_residue_select_all_policy(),
										align::common_atom_select_cb_policy(),
										post_proc,
										ssap_score_accuracy::HIGH,
										num,
										dist_score
									);
								}
							}
//						}
					}
//				}
//			}
		}
	}

	/// \brief Append all varieties of score object of the specified type to the specified ptr_vector of scores
	template <typename T>
	void score_variety_factory::leng_choice_and_subs_matrices_append(boost::ptr_vector<aligned_pair_score> &prm_scores ///< TODOCUMENT
	                                                                 ) {
		const boost::ptr_vector<length_getter> all_length_getters        = cath::score::get_all_length_getters();
		const substitution_matrix_vec          all_substitution_matrices = get_all_substitution_matrices();
		for (const length_getter &the_length_getter : all_length_getters) {
			for (const substitution_matrix &the_substitution_matrix : all_substitution_matrices) {
				boost::assign::ptr_push_back< T >( prm_scores )( the_length_getter, the_substitution_matrix );
			}
		}
	}

	/// \brief Append all varieties of score object of the specified type to the specified ptr_vector of scores
	template <typename T>
	void score_variety_factory::subs_matrices_append(boost::ptr_vector<aligned_pair_score> &prm_scores, ///< TODOCUMENT
	                                                 const length_getter                   &prm_length_getter ///< TODOCUMENT
	                                                 ) {
		const substitution_matrix_vec all_substitution_matrices = get_all_substitution_matrices();
		for (const substitution_matrix &the_substitution_matrix : all_substitution_matrices) {
			boost::assign::ptr_push_back< T >( prm_scores )( prm_length_getter, the_substitution_matrix );
		}
	}

	/// \brief TODOCUMENT
	///
	/// Default is to append a single, default-constructed T
	template <typename T>
	void score_variety_factory::append_all_varieties(boost::ptr_vector<aligned_pair_score> &prm_scores
	                                                 ) {
		boost::assign::ptr_push_back< T >( prm_scores )( );
	}

	template <>
	void score_variety_factory::append_all_varieties< drmsd_score               >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< gsas_score                >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< length_score              >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< lddt_score                >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< mi_score                  >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< rmsd_score                >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< sas_score                 >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< sequence_similarity_score >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< si_score                  >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< ssap_score                >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< structal_score            >( boost::ptr_vector<aligned_pair_score> & );
	template <>
	void score_variety_factory::append_all_varieties< tm_score                  >( boost::ptr_vector<aligned_pair_score> & );

	/// \brief TODOCUMENT
	template <typename T>
	boost::ptr_vector<aligned_pair_score> make_all_varieties() {
		boost::ptr_vector<aligned_pair_score> new_scores;
		score_variety_factory::append_all_varieties<T>(new_scores);
		return new_scores;
	}

} // namespace cath::score::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_DETAIL_ALIGNED_PAIR_SCORE_LIST_APPEND_HPP
