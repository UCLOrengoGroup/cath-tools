/// \file
/// \brief The true_false_pos_neg_list class definitions

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

#include "true_false_pos_neg_list.hpp"

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/score/true_pos_false_neg/classn_rate_stat.hpp"
#include "cath/score/true_pos_false_neg/classn_stat_pair_series.hpp"
#include "cath/score/true_pos_false_neg/true_false_pos_neg.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace std;

/// \brief Ctor from a vector of true_false_pos_neg objects
true_false_pos_neg_list::true_false_pos_neg_list(true_false_pos_neg_vec prm_tfpns ///< The vector of true_false_pos_neg objects with which this true_false_pos_neg_list should be constructed
                                                 ) : tfpns { std::move( prm_tfpns ) } {
}

/// \brief Return whether the vector is empty
bool true_false_pos_neg_list::empty() const {
	return tfpns.empty();
}

/// \brief Return the number of true_false_pos_neg objects
size_t true_false_pos_neg_list::size() const {
	return tfpns.size();
}

/// \brief Push another true_false_pos_neg object to the back of the list
void true_false_pos_neg_list::push_back(const true_false_pos_neg &prm_tfpn ///< The true_false_pos_neg object to push to the back of the true_false_pos_neg_list
                                        ) {
	tfpns.push_back( prm_tfpn );
}

/// \brief Subscript operator returning a const reference to the true_false_pos_neg object at the specified index
const true_false_pos_neg & true_false_pos_neg_list::operator[](const size_t &prm_index ///< The index of the true_false_pos_neg to return
                                                               ) const {
	return tfpns[ prm_index ];
}

/// \brief Return a begin() const_iterator
///
/// This is part of making the class into a (readonly) range
true_false_pos_neg_list::const_iterator true_false_pos_neg_list::begin() const {
	return tfpns.begin();
}

/// \brief Return an end() const_iterator
///
/// This is part of making the class into a (readonly) range
true_false_pos_neg_list::const_iterator true_false_pos_neg_list::end() const {
	return tfpns.end();
}

/// \brief TODOCUMENT
///
/// \relates true_false_pos_neg_list
size_rational cath::score::get_classn_stat_val_of_index(const true_false_pos_neg_list &prm_true_false_pos_neg_list, ///< TODOCUMENT
                                                        const classn_stat             &prm_classn_stat,              ///< TODOCUMENT
                                                        const size_t                  &prm_index                    ///< TODOCUMENT
                                                        ) {
	return prm_classn_stat.calculate( prm_true_false_pos_neg_list[ prm_index ] );
}

/// \brief TODOCUMENT
///
/// \relates true_false_pos_neg_list
size_rational_vec cath::score::get_classn_stat_vals(const true_false_pos_neg_list &prm_tfpns,      ///< TODOCUMENT
                                                    const classn_stat             &prm_classn_stat ///< TODOCUMENT
                                                    ) {
	return transform_build<size_rational_vec>(
		prm_tfpns,
		[&] (const true_false_pos_neg &x) { return prm_classn_stat.calculate( x ); }
	);
}

/// \brief TODOCUMENT
///
/// \relates true_false_pos_neg_list
classn_stat_pair_series cath::score::get_classn_stat_pair_series(const true_false_pos_neg_list &prm_tfpns,         ///< TODOCUMENT
                                                                 const string                  &prm_name,          ///< TODOCUMENT
                                                                 const classn_stat             &prm_classn_stat_a, ///< TODOCUMENT
                                                                 const classn_stat             &prm_classn_stat_b  ///< TODOCUMENT
                                                                 ) {
	return classn_stat_pair_series(
		transform_build<doub_doub_pair_vec>(
			prm_tfpns,
			[&] (const true_false_pos_neg &x) {
				return make_pair(
					calculate_and_convert( prm_classn_stat_a, x ),
					calculate_and_convert( prm_classn_stat_b, x )
				);
			}
		),
		prm_name
	);
}

/// \brief TODOCUMENT
///
/// \relates true_false_pos_neg_list
double cath::score::area_under_curve(const true_false_pos_neg_list &prm_tfpns,         ///< TODOCUMENT
                                     const classn_stat             &prm_classn_stat_x, ///< TODOCUMENT
                                     const classn_stat             &prm_classn_stat_y  ///< TODOCUMENT
                                     ) {
	return area_under_curve( get_classn_stat_pair_series( prm_tfpns, "", prm_classn_stat_x, prm_classn_stat_y ) );
}

/// \brief TODOCUMENT
///
/// \relates true_false_pos_neg_list
classn_stat_pair_series cath::score::get_roc_series(const true_false_pos_neg_list &prm_tfpns, ///< TODOCUMENT
                                                    const string                  &prm_name   ///< TODOCUMENT
                                                    ) {
	return get_classn_stat_pair_series(
		prm_tfpns, prm_name, roc_rates::first_type(), roc_rates::second_type()
	);
}

/// \brief TODOCUMENT
///
/// \relates true_false_pos_neg_list
double cath::score::area_under_roc_curve(const true_false_pos_neg_list &prm_tfpns ///< TODOCUMENT
                                         ) {
	return area_under_curve( get_roc_series( prm_tfpns, "" ) );
}
