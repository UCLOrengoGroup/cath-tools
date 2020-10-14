/// \file
/// \brief The classn_stat_pair_series_list class definitions

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

#include "classn_stat_pair_series_list.hpp"

#include <boost/range/algorithm/find_if.hpp>

#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/score/true_pos_false_neg/classn_stat_pair_series.hpp"

using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::range::find_if;

/// \brief TODOCUMENT
classn_stat_pair_series_list::classn_stat_pair_series_list(classn_stat_pair_series_vec prm_classn_stat_pair_serieses ///< TODOCUMENT
                                                           ) : classn_stat_pair_serieses{ std::move( prm_classn_stat_pair_serieses ) } {
}

/// \brief TODOCUMENT
bool classn_stat_pair_series_list::empty() const {
	return classn_stat_pair_serieses.empty();
}

/// \brief TODOCUMENT
size_t classn_stat_pair_series_list::size() const {
	return classn_stat_pair_serieses.size();
}

/// \brief TODOCUMENT
const classn_stat_pair_series & classn_stat_pair_series_list::operator[](const size_t &prm_index ///< TODOCUMENT
                                                                         ) const {
	return classn_stat_pair_serieses[ prm_index ];
}

/// \brief TODOCUMENT
classn_stat_pair_series_list::const_iterator classn_stat_pair_series_list::begin() const {
	return cath::common::cbegin( classn_stat_pair_serieses );
}

/// \brief TODOCUMENT
classn_stat_pair_series_list::const_iterator classn_stat_pair_series_list::end() const {
	return cath::common::cend( classn_stat_pair_serieses );
}

/// \brief TODOCUMENT
///
/// \relates classn_stat_pair_series_list
const classn_stat_pair_series & cath::score::classn_stat_pair_series_list_of_name(const classn_stat_pair_series_list &prm_classn_stat_pair_serieses, ///< TODOCUMENT
                                                                                  const string                       &prm_name                       ///< TODOCUMENT
                                                                                  ) {
	const auto find_itr = find_if(
		prm_classn_stat_pair_serieses,
		[&] (const classn_stat_pair_series &x) {
			return ( x.get_name() == prm_name );
		}
	);
	if ( find_itr == cath::common::cend( prm_classn_stat_pair_serieses ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find requested name \"" + prm_name + "\" in prm_classn_stat_pair_serieses"));
	}
	return *find_itr;
}
