/// \file
/// \brief The common_residue_select_min_score_policy class definitions

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

#include "common_residue_select_min_score_policy.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;

constexpr double common_residue_select_min_score_policy::MIN_CUTOFF;
constexpr double common_residue_select_min_score_policy::MAX_CUTOFF;
constexpr double common_residue_select_min_score_policy::DEFAULT_CUTOFF;

/// \brief TODOCUMENT
double common_residue_select_min_score_policy::get_score_cutoff() const {
	return score_cutoff;
}

/// \brief TODOCUMENT
size_vec common_residue_select_min_score_policy::do_select_common_residues_with_scores(const doub_doub_pair_vec &prm_scores ///< TODOCUMENT
                                                                                       ) const {
	size_vec the_indices;
	the_indices.reserve( prm_scores.size() );
	for (const size_t &index_ctr : indices( prm_scores.size() ) ) {
		const doub_doub_pair &score_pair = prm_scores[index_ctr];
		if ( min( score_pair.first, score_pair.second ) > get_score_cutoff() ) {
			the_indices.push_back( index_ctr );
		}
	}
	return the_indices;
}

/// \brief TODOCUMENT
string common_residue_select_min_score_policy::do_get_descriptive_name() const {
	return "select_min_score[" + lexical_cast<string>(get_score_cutoff()) + "]";
}

/// \brief TODOCUMENT
unique_ptr<common_residue_selection_policy> common_residue_select_min_score_policy::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Ctor for common_residue_select_min_score_policy
common_residue_select_min_score_policy::common_residue_select_min_score_policy(const double &prm_score_cutoff ///< TODOCUMENT
                                                                               ) : score_cutoff(prm_score_cutoff) {
	using boost::math::isfinite;
	if (!isfinite(get_score_cutoff())) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Argument score_cutoff must be a normal, finite floating-point number"));
	}
	if (get_score_cutoff() < MIN_CUTOFF || get_score_cutoff() > MAX_CUTOFF) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Argument score_cutoff is " + lexical_cast<string>(get_score_cutoff()) +
			" but should be between " + lexical_cast<string>(MIN_CUTOFF) +
			" and " + lexical_cast<string>(MAX_CUTOFF) +
			" (inclusive)"
		));
	}
}

/// \brief TODOCUMENT
bool common_residue_select_min_score_policy::do_less_than_with_same_dynamic_type(const common_residue_selection_policy &/*prm_common_residue_selection_policy*/ ///< TODOCUMENT
                                                                                 ) const {
	return false;
}

