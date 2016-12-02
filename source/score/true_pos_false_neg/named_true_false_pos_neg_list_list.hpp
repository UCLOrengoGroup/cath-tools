/// \file
/// \brief The named_true_false_pos_neg_list_list class header

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

#ifndef _CATH_TOOLS_SOURCE_SCORE_TRUE_POS_FALSE_NEG_NAMED_TRUE_FALSE_POS_NEG_LIST_LIST_H
#define _CATH_TOOLS_SOURCE_SCORE_TRUE_POS_FALSE_NEG_NAMED_TRUE_FALSE_POS_NEG_LIST_LIST_H

#include "common/type_aliases.hpp"
#include "score/score_type_aliases.hpp"

namespace cath { namespace score { class classn_stat; } }
namespace cath { namespace score { class classn_stat_pair_series_list; } }

namespace cath {
	namespace score {

		/// \brief TODOCUMENT
		class named_true_false_pos_neg_list_list final {
		private:
			/// \brief TODOCUMENT
			named_true_false_pos_neg_list_vec named_true_false_pos_neg_lists;

		public:
			using const_iterator = named_true_false_pos_neg_list_vec_citr;

			named_true_false_pos_neg_list_list(const named_true_false_pos_neg_list_vec &);

			bool empty() const;
			size_t size() const;

			const named_true_false_pos_neg_list & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		classn_stat_pair_series_list make_classn_stat_pair_series_list(const named_true_false_pos_neg_list_list &,
		                                                               const classn_stat &,
		                                                               const classn_stat &);

		str_doub_pair_vec areas_under_roc_curves(const named_true_false_pos_neg_list_list &);

	} // namespace score
} // namespace cath

#endif
