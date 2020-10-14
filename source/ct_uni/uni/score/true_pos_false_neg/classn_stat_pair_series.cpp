/// \file
/// \brief The classn_stat_pair_series class definitions

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

#include "classn_stat_pair_series.hpp"

#include <boost/range/numeric.hpp>

#include "common/algorithm/adjacent_accumulate.hpp"
#include "common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "score/true_pos_false_neg/classn_stat.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace std;

/// \brief TODOCUMENT
classn_stat_pair_series::classn_stat_pair_series(doub_doub_pair_vec prm_data, ///< TODOCUMENT
                                                 string             prm_name  ///< TODOCUMENT
                                                 ) : data { std::move( prm_data ) },
                                                     name { std::move( prm_name ) } {
}

/// \brief TODOCUMENT
bool classn_stat_pair_series::empty() const {
	return data.empty();
}

/// \brief TODOCUMENT
size_t classn_stat_pair_series::size() const {
	return data.size();
}

/// \brief TODOCUMENT
const doub_doub_pair & classn_stat_pair_series::operator[](const size_t &prm_index ///< TODOCUMENT
                                                           ) const {
	return data[ prm_index ];
}

/// \brief TODOCUMENT
classn_stat_pair_series::const_iterator classn_stat_pair_series::begin() const {
	return common::cbegin( data );
}

/// \brief TODOCUMENT
classn_stat_pair_series::const_iterator classn_stat_pair_series::end() const {
	return common::cend( data );
}

///// \brief TODOCUMENT
//const doub_doub_pair_vec & classn_stat_pair_series::get_doub_doub_data() const {
//	return data;
//}

/// \brief TODOCUMENT
const string & classn_stat_pair_series::get_name() const {
	return name;
}

/// \brief TODOCUMENT
///
/// \todo Write an "adjacented" adaptor that returns pairs of references to consecutive pairs of
///       elements in a range. Then use that here in tandem with the standard accumulate.
///
/// \relates classn_stat_pair_series
double cath::score::area_under_curve(const classn_stat_pair_series &prm_curve ///< TODOCUMENT
                                     ) {
	using value_type = range_value_t<const classn_stat_pair_series>;

	return adjacent_accumulate(
		prm_curve,
		0.0,
		[] (const value_type &a, const value_type &b) {
			// Sanity check the inputs...
			//
			// Check that the second x-axis value isn't greater
			// (but allow it to be equal because curves like ROC curves should be allowed to go straight up)
			if ( a.first > b.first ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot calculate area under curve that isn't monotonically increasing over the x-axis"));
			}
			// Check that neither y-axis value is less than 0
			if ( a.second < 0.0 || b.second < 0.0 ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot calculate area under curve that goes below the x-axis"));
			}

			// Calculate and return the area of the trapezium defined by dropping the two points down to the x-axis
			const auto x_diff =   b.first  - a.first;
			const auto y_mean = ( a.second + b.second ) / 2.0;
			return ( x_diff * y_mean );
		}
	);
}
