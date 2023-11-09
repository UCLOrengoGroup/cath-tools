/// \file
/// \brief The common_residue_select_best_score_percent_policy class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY_COMMON_RESIDUE_SELECT_BEST_SCORE_PERCENT_POLICY_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY_COMMON_RESIDUE_SELECT_BEST_SCORE_PERCENT_POLICY_HPP

#include "cath/alignment/common_residue_selection_policy/common_residue_score_based_selection_policy.hpp"

namespace cath::align {

	/// \brief TODOCUMENT
	class common_residue_select_best_score_percent_policy final : public common_residue_score_based_selection_policy {
	private:
		/// \brief TODOCUMENT
		double best_score_percentage;

		[[nodiscard]] size_vec    do_select_common_residues_with_scores( const doub_doub_pair_vec & ) const final;
		[[nodiscard]] std::string do_get_descriptive_name() const final;
		[[nodiscard]] std::unique_ptr<common_residue_selection_policy> do_clone() const final;

		[[nodiscard]] bool do_less_than_with_same_dynamic_type( const common_residue_selection_policy & ) const final;

	  public:
		explicit common_residue_select_best_score_percent_policy(const double & = DEFAULT_BEST_SCORE_PERCENTAGE);

		/// \brief TODOCUMENT
		static constexpr double MIN_BEST_SCORE_PERCENTAGE     =   0.0;

		/// \brief TODOCUMENT
		static constexpr double MAX_BEST_SCORE_PERCENTAGE     = 100.0;

		/// \brief TODOCUMENT
		static constexpr double DEFAULT_BEST_SCORE_PERCENTAGE =  70.0;
	};

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY_COMMON_RESIDUE_SELECT_BEST_SCORE_PERCENT_POLICY_HPP
