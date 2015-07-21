/// \file
/// \brief The common_residue_select_best_score_percent_policy class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "common_residue_select_best_score_percent_policy.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/stable_sort.hpp>
#include <boost/range.hpp>

#include "common/algorithm/copy_build.h"
#include "common/algorithm/sort_uniq_build.h"
#include "common/algorithm/transform_build.h"
#include "common/clone/make_uptr_clone.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"

#include <algorithm>
#include <numeric>

using namespace boost::math;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;
using boost::accumulate;
using boost::irange;
using boost::range::find_if;
using boost::range::sort;
using boost::range::stable_sort;
using boost::sub_range;

/// \brief TODOCUMENT
const double common_residue_select_best_score_percent_policy::MIN_BEST_SCORE_PERCENTAGE     = {   0.0 };

/// \brief TODOCUMENT
const double common_residue_select_best_score_percent_policy::MAX_BEST_SCORE_PERCENTAGE     = { 100.0 };

/// \brief TODOCUMENT
const double common_residue_select_best_score_percent_policy::DEFAULT_BEST_SCORE_PERCENTAGE = {  70.0 };

/// \brief TODOCUMENT
size_vec common_residue_select_best_score_percent_policy::do_select_common_residues_with_scores(const doub_doub_pair_vec &arg_scores ///< TODOCUMENT
                                                                                                ) const {
	// Grab the a list of the smaller value from each pair of scores
	const auto min_scores = transform_build<doub_vec>(
		arg_scores,
		[] (const doub_doub_pair &x) { return min( x.first, x.second ); }
	);

	// Build a stable_sorted list of indices in descending order of the score to which each corresponds
	const auto score_sorted_indices = stable_sort_build<size_vec>(
		irange( 0_z, min_scores.size() ),
		[&] (const size_t &x, const size_t &y) {
			return ( min_scores.at( x ) > min_scores.at( y ) );
		}
	);

	// Calculate the total score of all values and the cutoff for best_score_percentage % of it
	const double total_score            = accumulate( min_scores, 0.0 );
	const double percentile_total_score = total_score * best_score_percentage / 100.0;

	// Step through score_sorted_indices and grab them until the sum of their scores exceeds
	// percentile_total_score, then return a sorted copy of all those indices
	double temp_sum( 0.0 );
	return sort_build<size_vec>(
		common::cbegin( score_sorted_indices ),
		find_if(
			score_sorted_indices,
			[&] (const size_t &x) {
				temp_sum += min_scores[ x ];
				return ( temp_sum > percentile_total_score );
			}
		)
	);
}

/// \brief TODOCUMENT
string common_residue_select_best_score_percent_policy::do_get_descriptive_name() const {
	return "select_best_score_percent[" + lexical_cast<string>( best_score_percentage ) + "]";
}

/// \brief TODOCUMENT
unique_ptr<common_residue_selection_policy> common_residue_select_best_score_percent_policy::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Ctor for common_residue_select_best_score_percent_policy
common_residue_select_best_score_percent_policy::common_residue_select_best_score_percent_policy(const double &arg_best_score_percentage ///< TODOCUMENT
                                                                                                 ) : best_score_percentage(arg_best_score_percentage) {
	using boost::math::isfinite;
	if ( ! isfinite( best_score_percentage ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Argument best_score_percentage must be a normal, finite floating-point number"));
	}
	if ( best_score_percentage < MIN_BEST_SCORE_PERCENTAGE || best_score_percentage > MAX_BEST_SCORE_PERCENTAGE ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Argument best_score_percentage is " + lexical_cast<string>( best_score_percentage     ) +
			" but should be between "            + lexical_cast<string>( MIN_BEST_SCORE_PERCENTAGE ) +
			" and "                              + lexical_cast<string>( MAX_BEST_SCORE_PERCENTAGE ) +
			" (inclusive)"
		));
	}
}

/// \brief TODOCUMENT
bool common_residue_select_best_score_percent_policy::do_less_than_with_same_dynamic_type(const common_residue_selection_policy &/*arg_common_residue_selection_policy*/ ///< TODOCUMENT
                                                                                          ) const {
	return false;
}
