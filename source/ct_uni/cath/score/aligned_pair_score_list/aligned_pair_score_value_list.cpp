/// \file
/// \brief The aligned_pair_score_value_list class definitions

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

#include "aligned_pair_score_value_list.hpp"

#include <iostream> // ***** TEMPORARY *****

#include <boost/lexical_cast.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/range/algorithm/adjacent_find.hpp>

#include <spdlog/spdlog.h>

#include "cath/alignment/alignment.hpp"
#include "cath/common/algorithm/sort_uniq_copy.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/boost_addenda/tribool/tribool.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/score/aligned_pair_score/aligned_pair_score.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::cath::score;
using namespace ::std;

using ::boost::lexical_cast;
using ::boost::property_tree::ptree;
using ::boost::range::adjacent_find;
using ::boost::tribool;

/// \brief TODOCUMENT
void aligned_pair_score_value_list::add_score_and_value(const aligned_pair_score &prm_score, ///< TODOCUMENT
                                                        const score_value        &prm_value  ///< TODOCUMENT
                                                        ) {
	score_values.push_back(prm_value);
	scores.add_aligned_pair_score(prm_score);
}

/// \brief TODOCUMENT
size_t aligned_pair_score_value_list::size() const {
	return score_values.size();
}

/// \brief TODOCUMENT
score_value aligned_pair_score_value_list::get_value_of_index(const size_t &prm_index ///< TODOCUMENT
                                                              ) const {
	return score_values[ prm_index ];
}

/// \brief TODOCUMENT
const aligned_pair_score & aligned_pair_score_value_list::get_aligned_pair_score_of_index(const size_t &prm_index ///< TODOCUMENT
                                                                                          ) const {
	return scores[ prm_index ];
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
void cath::score::warn_on_duplicate_human_friendly_names(const aligned_pair_score_value_list &prm_aligned_pair_score_value_list ///< TODOCUMENT
                                                         ) {
	const size_t num_scores = prm_aligned_pair_score_value_list.size();
	const auto sorted_human_friendly_names = sort_copy( transform_build<str_vec>(
		indices( num_scores ),
		[&] (const size_t &x) {
			const auto &score = prm_aligned_pair_score_value_list.get_aligned_pair_score_of_index( x );
			return score.human_friendly_short_name();
		}
	) );
	const auto adjacent_itr = adjacent_find( sorted_human_friendly_names );
	if ( adjacent_itr != cend( sorted_human_friendly_names ) ) {
		::spdlog::warn( "aligned_pair_score_value_list contains duplicates (eg {})", *adjacent_itr );
	}
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
str_vec cath::score::get_all_names(const aligned_pair_score_value_list &prm_aligned_pair_score_value_list ///< TODOCUMENT
                                   ) {
	const size_t num_score_values = prm_aligned_pair_score_value_list.size();
	return transform_build<str_vec>(
		indices( num_score_values ),
		[&] (const size_t &x) { return prm_aligned_pair_score_value_list.get_aligned_pair_score_of_index( x ).human_friendly_short_name(); }
	);
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
void cath::score::save_to_ptree(ptree                               &prm_ptree,           ///< TODOCUMENT
                                const aligned_pair_score_value_list &prm_score_value_list ///< TODOCUMENT
                                ) {
	const size_t num_scores = prm_score_value_list.size();
	prm_ptree.put_child( "scores", ptree() );
	auto &ptree_scores_array = prm_ptree.get_child( "scores" );

	for (const size_t &score_ctr : indices( num_scores ) ) {
		const score_value         the_score_value           = prm_score_value_list.get_value_of_index( score_ctr );
		const aligned_pair_score &the_score_type            = prm_score_value_list.get_aligned_pair_score_of_index( score_ctr );
		const string              human_friendly_short_name = the_score_type.human_friendly_short_name();
		const tribool             higher_is_better          = the_score_type.higher_is_better();

//		if ( ! indeterminate( higher_is_better ) ) {
			ptree score_ptree;

			score_ptree.put( "score_type",       human_friendly_short_name );
			score_ptree.put( "score_value",      the_score_value           );
			ostringstream higher_is_better_ss;
			higher_is_better_ss << boolalpha << higher_is_better;
			score_ptree.put( "higher_is_better", higher_is_better_ss.str() );

//			cerr << "Processing score : "    << human_friendly_short_name
//			     << " with value "           << the_score_value
//			     << " and higher_is_better " << higher_is_better_ss.str()
//			     << endl;

			// Previously just used `prm_ptree.add_child("scores.", score_ptree );` but this caused the debug version to fail
			// failed an assert with message: `Empty path not allowed for put_child.`
			// See: http://stackoverflow.com/a/8187487
			ptree_scores_array.push_back( make_pair( "", score_ptree ) );

		}
//	}
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
aligned_pair_score_value_list cath::score::make_aligned_pair_score_value_list(const aligned_pair_score_list &prm_scores,    ///< TODOCUMENT
                                                                              const alignment               &prm_alignment, ///< TODOCUMENT
                                                                              const protein                 &prm_protein_a, ///< TODOCUMENT
                                                                              const protein                 &prm_protein_b  ///< TODOCUMENT
                                                                              ) {
	// Sanity-check that there are two entries in the alignment
	if ( 2 != prm_alignment.num_entries() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to make_aligned_pair_score_value_list() for alignment that contains "
			+ lexical_cast<string>( prm_alignment.num_entries() )
			+ " but should contain 2"
		));
	}

	// Create a new aligned_pair_score_value_list
	aligned_pair_score_value_list new_aligned_pair_score_value_list;

	// Loop over the scores
	const size_t num_scores = prm_scores.size();
	for (const size_t &score_ctr : indices( num_scores ) ) {
		const aligned_pair_score &score = prm_scores[score_ctr];

		// Calculate the value for the score and append it to the new aligned_pair_score_value_list
		const score_value value = score.calculate(prm_alignment, prm_protein_a, prm_protein_b);
		new_aligned_pair_score_value_list.add_score_and_value( score, value );
	}

	// Return the newly created aligned_pair_score_value_list
	return new_aligned_pair_score_value_list;
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
ostream & cath::score::operator<<(ostream                             &prm_os,                           ///< TODOCUMENT
                                  const aligned_pair_score_value_list &prm_aligned_pair_score_value_list ///< TODOCUMENT
                                  ) {
	const size_t num_scores = prm_aligned_pair_score_value_list.size();
	for (const size_t &score_ctr : indices( num_scores ) ) {
		const aligned_pair_score &score = prm_aligned_pair_score_value_list.get_aligned_pair_score_of_index(score_ctr);
		const score_value         value = prm_aligned_pair_score_value_list.get_value_of_index(score_ctr);
		prm_os << score.human_friendly_short_name();
		prm_os << ":";
		prm_os << value;
		prm_os << ";";
		prm_os << ( is_not_false( score.higher_is_better() ) ? value : -value );
		prm_os << "\n";
	}
	return prm_os;
}
