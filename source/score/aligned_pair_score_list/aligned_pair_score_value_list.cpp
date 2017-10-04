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

#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/range/algorithm/adjacent_find.hpp>

#include "alignment/alignment.hpp"
#include "common/algorithm/sort_uniq_copy.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/boost_addenda/tribool/tribool.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "score/aligned_pair_score/aligned_pair_score.hpp"

#include <iostream> // ***** TEMPORARY *****

using namespace boost::logic;
using namespace boost::property_tree;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::lexical_cast;
using boost::range::adjacent_find;

/// \brief TODOCUMENT
void aligned_pair_score_value_list::add_score_and_value(const aligned_pair_score &arg_score, ///< TODOCUMENT
                                                        const score_value        &arg_value  ///< TODOCUMENT
                                                        ) {
	score_values.push_back(arg_value);
	scores.add_aligned_pair_score(arg_score);
}

/// \brief TODOCUMENT
size_t aligned_pair_score_value_list::size() const {
	return score_values.size();
}

/// \brief TODOCUMENT
score_value aligned_pair_score_value_list::get_value_of_index(const size_t &arg_index ///< TODOCUMENT
                                                              ) const {
	return score_values[ arg_index ];
}

/// \brief TODOCUMENT
const aligned_pair_score & aligned_pair_score_value_list::get_aligned_pair_score_of_index(const size_t &arg_index ///< TODOCUMENT
                                                                                          ) const {
	return scores[ arg_index ];
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
void cath::score::warn_on_duplicate_human_friendly_names(const aligned_pair_score_value_list &arg_aligned_pair_score_value_list ///< TODOCUMENT
                                                         ) {
	const size_t num_scores = arg_aligned_pair_score_value_list.size();
	const auto sorted_human_friendly_names = sort_copy( transform_build<str_vec>(
		indices( num_scores ),
		[&] (const size_t &x) {
			const auto &score = arg_aligned_pair_score_value_list.get_aligned_pair_score_of_index( x );
			return score.human_friendly_short_name();
		}
	) );
	const auto adjacent_itr = adjacent_find( sorted_human_friendly_names );
	if ( adjacent_itr != common::cend( sorted_human_friendly_names ) ) {
		BOOST_LOG_TRIVIAL( warning ) << "aligned_pair_score_value_list contains duplicates (eg " << *adjacent_itr << ")";
	}
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
str_vec cath::score::get_all_names(const aligned_pair_score_value_list &arg_aligned_pair_score_value_list ///< TODOCUMENT
                                   ) {
	const size_t num_score_values = arg_aligned_pair_score_value_list.size();
	return transform_build<str_vec>(
		indices( num_score_values ),
		[&] (const size_t &x) { return arg_aligned_pair_score_value_list.get_aligned_pair_score_of_index( x ).human_friendly_short_name(); }
	);
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
void cath::score::save_to_ptree(ptree                               &arg_ptree,           ///< TODOCUMENT
                                const aligned_pair_score_value_list &arg_score_value_list ///< TODOCUMENT
                                ) {
	const size_t num_scores = arg_score_value_list.size();
	arg_ptree.put_child( "scores", ptree() );
	auto &ptree_scores_array = arg_ptree.get_child( "scores" );

	for (size_t score_ctr = 0; score_ctr < num_scores; ++score_ctr) {
		arg_score_value_list.get_value_of_index(score_ctr);
		const score_value         the_score_value           = arg_score_value_list.get_value_of_index( score_ctr );
		const aligned_pair_score &the_score_type            = arg_score_value_list.get_aligned_pair_score_of_index( score_ctr );
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

			// Previously just used `arg_ptree.add_child("scores.", score_ptree );` but this caused the debug version to fail
			// failed an assert with message: `Empty path not allowed for put_child.`
			// See: http://stackoverflow.com/a/8187487
			ptree_scores_array.push_back( make_pair( "", score_ptree ) );

		}
//	}
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
aligned_pair_score_value_list cath::score::make_aligned_pair_score_value_list(const aligned_pair_score_list &arg_scores,    ///< TODOCUMENT
                                                                              const alignment               &arg_alignment, ///< TODOCUMENT
                                                                              const protein                 &arg_protein_a, ///< TODOCUMENT
                                                                              const protein                 &arg_protein_b  ///< TODOCUMENT
                                                                              ) {
	// Sanity-check that there are two entries in the alignment
	if ( 2 != arg_alignment.num_entries() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to make_aligned_pair_score_value_list() for alignment that contains "
			+ lexical_cast<string>( arg_alignment.num_entries() )
			+ " but should contain 2"
		));
	}

	// Create a new aligned_pair_score_value_list
	aligned_pair_score_value_list new_aligned_pair_score_value_list;

	// Loop over the scores
	const size_t num_scores = arg_scores.size();
	for (size_t score_ctr = 0; score_ctr < num_scores; ++score_ctr) {
		const aligned_pair_score &score = arg_scores[score_ctr];

		// Calculate the value for the score and append it to the new aligned_pair_score_value_list
		const score_value value = score.calculate(arg_alignment, arg_protein_a, arg_protein_b);
		new_aligned_pair_score_value_list.add_score_and_value( score, value );
	}

	// Return the newly created aligned_pair_score_value_list
	return new_aligned_pair_score_value_list;
}

/// \brief TODOCUMENT
///
/// \relates aligned_pair_score_value_list
ostream & cath::score::operator<<(ostream                             &arg_os,                           ///< TODOCUMENT
                                  const aligned_pair_score_value_list &arg_aligned_pair_score_value_list ///< TODOCUMENT
                                  ) {
	const size_t num_scores = arg_aligned_pair_score_value_list.size();
	for (size_t score_ctr = 0; score_ctr < num_scores; ++score_ctr) {
		const aligned_pair_score &score = arg_aligned_pair_score_value_list.get_aligned_pair_score_of_index(score_ctr);
		const score_value         value = arg_aligned_pair_score_value_list.get_value_of_index(score_ctr);
		arg_os << score.human_friendly_short_name();
		arg_os << ":";
		arg_os << value;
		arg_os << ";";
		arg_os << ( is_not_false( score.higher_is_better() ) ? value : -value );
		arg_os << "\n";
	}
	return arg_os;
}
