/// \file
/// \brief The true_false_pos_neg_list class header

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

#ifndef _CATH_TOOLS_SOURCE_SCORE_TRUE_POS_FALSE_NEG_TRUE_FALSE_POS_NEG_LIST_H
#define _CATH_TOOLS_SOURCE_SCORE_TRUE_POS_FALSE_NEG_TRUE_FALSE_POS_NEG_LIST_H

#include "score/score_type_aliases.hpp"
#include "score/true_pos_false_neg/classn_stat.hpp"

namespace cath { namespace score { class classn_stat_pair_series; } }

namespace cath {
	namespace score {

		/// \brief A list of true_false_pos_neg entries, typically associated with the
		///        results
		///
		/// A ROC curve can be directly derived from a true_false_pos_neg_list but
		/// that involves a loss of information (because the true_false_pos_neg_list
		/// contains all absolute TP/FP/TN/FN counts)
		class true_false_pos_neg_list final {
		private:
			/// \brief The list of true_false_pos_neg objects
			true_false_pos_neg_vec tfpns;

		public:
			using const_iterator = true_false_pos_neg_vec_citr;

			explicit true_false_pos_neg_list(true_false_pos_neg_vec);

			bool empty() const;
			size_t size() const;

			void push_back(const true_false_pos_neg &);

			const true_false_pos_neg & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		size_rational get_classn_stat_val_of_index(const true_false_pos_neg_list &,
		                                           const classn_stat &,
		                                           const size_t &);

		size_rational_vec get_classn_stat_vals(const true_false_pos_neg_list &,
		                                       const classn_stat &);

		classn_stat_pair_series get_classn_stat_pair_series(const true_false_pos_neg_list &,
		                                                    const std::string &,
		                                                    const classn_stat &,
		                                                    const classn_stat &);

		double area_under_curve(const true_false_pos_neg_list &,
		                        const classn_stat &,
		                        const classn_stat &);

		classn_stat_pair_series get_roc_series(const true_false_pos_neg_list &,
		                                       const std::string &);

		double area_under_roc_curve(const true_false_pos_neg_list &);

	} // namespace score
} // namespace cath

#endif
