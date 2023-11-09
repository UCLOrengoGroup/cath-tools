/// \file
/// \brief The named_true_false_pos_neg_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_NAMED_TRUE_FALSE_POS_NEG_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_NAMED_TRUE_FALSE_POS_NEG_LIST_HPP

#include "cath/score/true_pos_false_neg/true_false_pos_neg_list.hpp"

namespace cath::score {

	/// \brief Wrap a true_false_pos_neg_list (ie a series of ROC-like data points) with
	///        a name (eg a series title for the algorithm/metric that generated the data)
	class named_true_false_pos_neg_list final {
	private:
		/// \brief The series of TFPN data points (from which a ROC-like curve can be derived)
		true_false_pos_neg_list the_list;

		/// \brief The name for the series of data
		std::string name;

	public:
		named_true_false_pos_neg_list(true_false_pos_neg_list,
		                              std::string);
		named_true_false_pos_neg_list(true_false_pos_neg_vec,
		                              std::string);

		[[nodiscard]] const true_false_pos_neg_list &get_list() const;
		[[nodiscard]] const std::string &            get_name() const;
	};

	size_rational get_classn_stat_val_of_index(const named_true_false_pos_neg_list &,
	                                           const classn_stat &,
	                                           const size_t &);

	size_rational_vec get_classn_stat_vals(const named_true_false_pos_neg_list &,
	                                       const classn_stat &);

	classn_stat_pair_series get_classn_stat_pair_series(const named_true_false_pos_neg_list &,
	                                                    const classn_stat &,
	                                                    const classn_stat &);

	double area_under_curve(const named_true_false_pos_neg_list &,
	                        const classn_stat &,
	                        const classn_stat &);

	classn_stat_pair_series get_roc_series(const named_true_false_pos_neg_list &);

	double area_under_roc_curve(const named_true_false_pos_neg_list &);

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_TRUE_POS_FALSE_NEG_NAMED_TRUE_FALSE_POS_NEG_LIST_HPP
