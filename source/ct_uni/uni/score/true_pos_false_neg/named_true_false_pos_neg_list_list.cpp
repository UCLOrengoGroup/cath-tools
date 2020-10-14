/// \file
/// \brief The named_true_false_pos_neg_list_list class definitions

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

#include "named_true_false_pos_neg_list_list.hpp"

#include "common/algorithm/transform_build.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "score/true_pos_false_neg/classn_stat_pair_series.hpp"
#include "score/true_pos_false_neg/classn_stat_pair_series_list.hpp"
#include "score/true_pos_false_neg/named_true_false_pos_neg_list.hpp"
#include "score/true_pos_false_neg/true_false_pos_neg.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::score;
// using namespace std;

/// \brief TODOCUMENT
named_true_false_pos_neg_list_list::named_true_false_pos_neg_list_list(named_true_false_pos_neg_list_vec prm_named_true_false_pos_neg_list_vec ///< TODOCUMENT
                                                                       ) : named_true_false_pos_neg_lists{ std::move( prm_named_true_false_pos_neg_list_vec ) } {

}

/// \brief TODOCUMENT
bool named_true_false_pos_neg_list_list::empty() const {
	return named_true_false_pos_neg_lists.empty();
}

/// \brief TODOCUMENT
size_t named_true_false_pos_neg_list_list::size() const {
	return named_true_false_pos_neg_lists.size();
}

/// \brief TODOCUMENT
const named_true_false_pos_neg_list & named_true_false_pos_neg_list_list::operator[](const size_t &prm_index ///< TODOCUMENT
                                                                                     ) const {
	return named_true_false_pos_neg_lists[ prm_index ];
}

/// \brief TODOCUMENT
named_true_false_pos_neg_list_list::const_iterator named_true_false_pos_neg_list_list::begin() const {
	return cath::common::cbegin( named_true_false_pos_neg_lists );
}

/// \brief TODOCUMENT
named_true_false_pos_neg_list_list::const_iterator named_true_false_pos_neg_list_list::end() const {
	return cath::common::cend( named_true_false_pos_neg_lists );
}

/// \brief TODOCUMENT
classn_stat_pair_series_list cath::score::make_classn_stat_pair_series_list(const named_true_false_pos_neg_list_list &prm_named_true_false_pos_neg_list_list, ///< TODOCUMENT
                                                                            const classn_stat                        &prm_classn_stat_a,                      ///< TODOCUMENT
                                                                            const classn_stat                        &prm_classn_stat_b                       ///< TODOCUMENT
                                                                            ) {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return classn_stat_pair_series_list{
		transform_build<classn_stat_pair_series_vec>(
			prm_named_true_false_pos_neg_list_list,
			[&] (const named_true_false_pos_neg_list &x) {
				return get_classn_stat_pair_series( x, prm_classn_stat_a, prm_classn_stat_b );
			}
		)
	};
}

/// \brief TODOCUMENT
str_doub_pair_vec cath::score::areas_under_roc_curves(const named_true_false_pos_neg_list_list &prm_named_true_false_pos_neg_list_list
                                                      ) {
	return transform_build<str_doub_pair_vec>(
		prm_named_true_false_pos_neg_list_list,
		[] (const named_true_false_pos_neg_list &x) {
			return make_pair( x.get_name(), area_under_roc_curve( x ) );
		}
	);
}

