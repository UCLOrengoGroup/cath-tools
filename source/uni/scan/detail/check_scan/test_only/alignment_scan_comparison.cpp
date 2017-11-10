/// \file
/// \brief The alignment_scan_comparison class definitions

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

#include "alignment_scan_comparison.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/map.hpp>

#include "common/algorithm/contains.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/sort_uniq_copy.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/type_aliases.hpp"

#include <iterator> // for end, inserter
#include <ostream>  // for ostream, operator<<, etc
#include <string>

using namespace cath;
using namespace cath::common;
using namespace cath::scan::detail;
using namespace std;

using boost::adaptors::map_keys;
using boost::algorithm::join;

/// \brief TODOCUMENT
alignment_scan_comparison & alignment_scan_comparison::operator+=(const pair<quad_and_rep_criteria_result, double> &arg_result_and_score ///< TODOCUMENT
                                                                  ) {
	const auto find_itr = score_by_result.find( arg_result_and_score.first );
	if ( find_itr == common::cend( score_by_result ) ) {
		score_by_result.insert( arg_result_and_score );
	}
	else {
		find_itr->second += arg_result_and_score.second;
	}
	return *this;
}

/// \brief TODOCUMENT
bool alignment_scan_comparison::has_score_of_criteria_result(const quad_and_rep_criteria_result &arg_comp ///< TODOCUMENT
                                                             ) const {
	return contains( score_by_result, arg_comp );
}

/// \brief TODOCUMENT
double alignment_scan_comparison::get_score_of_criteria_result(const quad_and_rep_criteria_result &arg_comp ///< TODOCUMENT
                                                               ) const {
	return score_by_result.at( arg_comp );
}

/// \brief TODOCUMENT
auto alignment_scan_comparison::begin() const -> const_iterator {
	return common::cbegin( score_by_result );
}

/// \brief TODOCUMENT
auto alignment_scan_comparison::end() const -> const_iterator {
	return common::cend  ( score_by_result );
}

/// \brief TODOCUMENT
///
/// \relates alignment_scan_comparison
vector<quad_and_rep_criteria_result> cath::scan::detail::get_criteria_results(const alignment_scan_comparison &arg_comp ///< TODOCUMENT
                                                                              ) {
	return copy_build<vector<quad_and_rep_criteria_result>>( arg_comp | map_keys );
}

/// \brief TODOCUMENT
///
/// \relates alignment_scan_comparison
vector<quad_and_rep_criteria_result> cath::scan::detail::get_criteria_results_sorted_by_score_desc(const alignment_scan_comparison &arg_comp ///< TODOCUMENT
                                                                                                   ) {
	return sort_copy(
		get_criteria_results( arg_comp ),
		[&] (const quad_and_rep_criteria_result &x, const quad_and_rep_criteria_result &y) {
			return arg_comp.get_score_of_criteria_result( x ) > arg_comp.get_score_of_criteria_result( y );
		}
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment_scan_comparison
ostream & cath::scan::detail::operator<<(ostream                         &arg_os,  ///< TODOCUMENT
                                         const alignment_scan_comparison &arg_comp ///< TODOCUMENT
                                         ) {
	const auto crit_results = get_criteria_results_sorted_by_score_desc( arg_comp );
	const auto the_strings = transform_build<str_vec>(
		crit_results,
		[&] (const quad_and_rep_criteria_result &crit_result) {
			return to_string( crit_result )
			       + " -> "
			       + std::to_string( arg_comp.get_score_of_criteria_result( crit_result ) );
		}
	);
	arg_os << (
		string( "alignment_scan_comparison[\n\t" )
		+
		join( the_strings, "\n\t" )
		+
		"\n]"
	);
	return arg_os;
}
