/// \file
/// \brief The aligned_pair_score_list class definitions

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

#include "aligned_pair_score_list.hpp"

#include <boost/log/trivial.hpp>
#include <boost/range/algorithm/adjacent_find.hpp>

#include "cath/common/algorithm/sort_uniq_copy.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/ptr_container/unique_ptr_functions.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/score/aligned_pair_score/aligned_pair_score.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::ptr_vector;
using boost::range::adjacent_find;

/// \brief Ctor for aligned_pair_score_list that allows the caller to specify a ptr_vector of aligned_pair_scores
aligned_pair_score_list::aligned_pair_score_list(const ptr_vector<aligned_pair_score> &prm_aligned_pair_scores ///< The ptr_vector of aligned_pair_scores with which this aligned_pair_score_list should be initialised
                                                 ) : scores(prm_aligned_pair_scores) {
}

/// \brief Add the specified aligned_pair_score to this list
void aligned_pair_score_list::add_aligned_pair_score(const aligned_pair_score &prm_aligned_pair_score ///< The aligned_pair_score to add to this list
                                                     ) {
	cath::common::push_back( scores, prm_aligned_pair_score.clone() );
}

/// \brief Return the number of scores currently held in this aligned_pair_score_list
///
/// This is a pass-through to size() method of scores
size_t aligned_pair_score_list::size() const {
	return scores.size();
}

/// \brief Return the score for the specified index
///
/// This is a pass-through to non-const subscript operator of scores
aligned_pair_score & aligned_pair_score_list::operator[](const size_t &prm_index ///< The index of the score to be retrieved
                                                         ) {
	return scores[prm_index];
}

/// \brief Return the score for the specified index
///
/// This is a pass-through to const subscript operator of scores
const aligned_pair_score & aligned_pair_score_list::operator[](const size_t &prm_index ///< The index of the score to be retrieved
                                                               ) const {
	return scores[prm_index];
}

/// \brief Standard begin() operator to allow aligned_pair_score_list to be accessed as range
aligned_pair_score_list::const_iterator aligned_pair_score_list::begin() const {
	return common::cbegin( scores );
}

/// \brief Standard end() operator to allow aligned_pair_score_list to be accessed as range
aligned_pair_score_list::const_iterator aligned_pair_score_list::end() const {
	return common::cend( scores );
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_list
void cath::score::warn_on_duplicate_human_friendly_names(const aligned_pair_score_list &prm_aligned_pair_score_list ///< TODOCUMENT
                                                         ) {
	const auto sorted_human_friendly_names = sort_copy( transform_build<str_vec>(
		prm_aligned_pair_score_list,
		[&] (const aligned_pair_score &x) {
			return x.human_friendly_short_name();
		}
	) );
	const auto adjacent_itr = adjacent_find( sorted_human_friendly_names );
	if ( adjacent_itr != common::cend( sorted_human_friendly_names ) ) {
		BOOST_LOG_TRIVIAL( warning ) << "aligned_pair_score_value_list contains duplicates (eg " << *adjacent_itr << ")";
	}
}

/// \brief Return a vector of the short names of all the aligned_pair_scores in the specified aligned_pair_score_list
///
/// \relates aligned_pair_score_list
str_vec cath::score::get_short_names(const aligned_pair_score_list &prm_aligned_pair_score_list ///< The aligned_pair_score_list from which the short names should be retrieved
                                     ) {
	return transform_build<str_vec>(
		prm_aligned_pair_score_list,
		[] (const aligned_pair_score &x) { return x.human_friendly_short_name(); }
	);
}

/// \brief Return a vector of the long names of all the aligned_pair_scores in the specified aligned_pair_score_list
///
/// \relates aligned_pair_score_list
str_vec cath::score::get_long_names(const aligned_pair_score_list &prm_aligned_pair_score_list ///< The aligned_pair_score_list from which the long names should be retrieved
                                    ) {
	return transform_build<str_vec>(
		prm_aligned_pair_score_list,
		[] (const aligned_pair_score &x) { return x.long_name(); }
	);
}

/// \brief Return a vector of the descriptions of all the aligned_pair_scores in the specified aligned_pair_score_list
///
/// \relates aligned_pair_score_list
str_vec cath::score::get_descriptions(const aligned_pair_score_list &prm_aligned_pair_score_list ///< The aligned_pair_score_list from which the descriptions should be retrieved
                                      ) {
	return transform_build<str_vec>(
		prm_aligned_pair_score_list,
		[] (const aligned_pair_score &x) { return x.description(); }
	);
}

/// \brief Return a vector of the references of all the aligned_pair_scores in the specified aligned_pair_score_list
///        (including empty strings for the aligned_pair_scores that have no references)
///
/// \relates aligned_pair_score_list
str_vec cath::score::get_references(const aligned_pair_score_list &prm_aligned_pair_score_list ///< The aligned_pair_score_list from which the references should be retrieved
		                            ) {
	return transform_build<str_vec>(
		prm_aligned_pair_score_list,
		[] (const aligned_pair_score &x) { return x.reference(); }
	);
}

// SSAP score:
//       pre_globalised__logged__ssap_score (current default)
//      post_globalised__logged__ssap_score
//  post_comparedalised__logged__ssap_score
//   pre_comparedalised__logged__ssap_score
//                local__logged__ssap_score
//         globalised__unlogged__ssap_score
//          middlised__unlogged__ssap_score
//              local__unlogged__ssap_score
// // Maximum number residues pairs in alignment that can be superposed with RMSD <= 3.0 A [requested by Christine in SCOP/CATH mapping workshop on 5th February 2013]
//
// GDT (global distance test)       [http://nar.oxfordjournals.org/content/31/13/3370.long http://en.wikipedia.org/wiki/Global_distance_test http://predictioncenter.org/casp/casp7/public/doc/LCS_GDT.README http://onlinelibrary.wiley.com/doi/10.1002/prot.21761/full]
// LCS (Longest Continuous Segment) [http://nar.oxfordjournals.org/content/31/13/3370.long http://en.wikipedia.org/wiki/Global_distance_test http://predictioncenter.org/casp/casp7/public/doc/LCS_GDT.README http://onlinelibrary.wiley.com/doi/10.1002/prot.21761/full]
//
// LCS_1_over_aligned
// LCS_1_over_longer
// LCS_1_over_shorter
// LCS_2_over_aligned
// LCS_2_over_longer
// LCS_2_over_shorter
// LCS_5_over_aligned
// LCS_5_over_longer
// LCS_5_over_shorter
// GDT_0.5_over_aligned
// GDT_0.5_over_longer
// GDT_0.5_over_shorter
// GDT_1.0_over_aligned
// GDT_1.0_over_longer
// GDT_1.0_over_shorter
// ...
// ...
// ...
// GDT_10.0_over_aligned
// GDT_10.0_over_longer
// GDT_10.0_over_shorter
// GDT_TS_XXX (global distance test total score,   an average of GDT_1.0_XXX, GDT_2.0_XXX, GDT_4.0_XXX and GDT_8.0_XXX)
// GDT_HA_XXX (global distance test high accuracy, an average of GDT_0.5_XXX, GDT_1.0_XXX, GDT_2.0_XXX and GDT_4.0_XXX)
//
// Quite a few others from Liisa Holm paper http://www.sciencedirect.com/science/article/pii/S0959440X09000621
//
// DRID (http://pubs.acs.org/doi/abs/10.1021/ct3003145) Complicated RMS between two vectors containing mean, root of variance and cube root of skew
// 		                                                of inverse distances to all atoms from some specific set of centroid atoms
//
// TM-Score
// sum_over_aligned_cas( 1 / (1 + d_i ^ 2 / d_0 ^ 2) ) / num_aligned_res
// (where d_i is distance between i_th pair of aligned CAs and d_0 := 1.24 * cube_root(num_aligned_res - 15) - 1.8
// (http://en.wikipedia.org/wiki/Template_modeling_score)
