/// \file
/// \brief The filter_vs_full_score class definitions

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

#include "filter_vs_full_score.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/type_aliases.hpp"

#include <iostream> // ***** TEMPORARY *****

using namespace boost::algorithm;
using namespace cath::common;
using namespace cath::index::filter;
using namespace cath::score;
using namespace std;

using boost::algorithm::is_any_of;
using boost::algorithm::starts_with;
using boost::algorithm::token_compress_on;

const string filter_vs_full_score::STRING_PREFIX( "filter_vs_full_score[" );
const string filter_vs_full_score::STRING_SUFFIX( "]"                     );

/// \brief Ctor for filter_vs_full_score
filter_vs_full_score::filter_vs_full_score(const double &arg_filter_score, ///< The filter score with which this should be constructed
                                           const double &arg_full_score    ///< The full score with which this should be constructed
                                           ) : filter_score( arg_filter_score ),
                                               full_score  ( arg_full_score   ) {
}

/// \brief Getter for the filter score
const double & filter_vs_full_score::get_filter_score() const {
	return filter_score;
}

/// \brief Getter for the full score
const double & filter_vs_full_score::get_full_score() const {
	return full_score;
}

/// \brief Standard equality operator for filter_vs_full_score
///
/// \relates filter_vs_full_score
bool cath::index::filter::operator==(const filter_vs_full_score &arg_filter_vs_full_score_a, ///< The first filter_vs_full_score to compare
                                     const filter_vs_full_score &arg_filter_vs_full_score_b  ///< The second filter_vs_full_score to compare
                                     ) {
	if ( arg_filter_vs_full_score_a.get_filter_score() != arg_filter_vs_full_score_b.get_filter_score() ) {
		return false;
	}
	if ( arg_filter_vs_full_score_a.get_full_score()   != arg_filter_vs_full_score_b.get_full_score()   ) {
		return false;
	}
	return true;
}

/// \brief Assess a real pair of filter/full scores against a filter attempt (which tries to identify
///        all full scores >= some value by selecting all filter scores >= some value)
///
/// \relates filter_vs_full_score
classn_outcome cath::index::filter::assess_real_scores_on_filter_attempt(const filter_vs_full_score &arg_real_scores,   ///< The real scores to be assessed
                                                                         const filter_vs_full_score &arg_filter_attempt ///< The filter attempt (interpretation: this attempts to identify all entries with full score >= filter_vs_full_score's full score by selecting all entries with filter score >= filter_vs_full_score's filter score)
                                                                         ) {
	const bool correct_answer     = ( arg_real_scores.get_full_score()   >= arg_filter_attempt.get_full_score()   );
	const bool decision_of_filter = ( arg_real_scores.get_filter_score() >= arg_filter_attempt.get_filter_score() );
	return outcome_of_correct_and_decision( correct_answer, decision_of_filter );
}

/// \brief Simple extraction operator for filter_vs_full_score
///
/// \relates filter_vs_full_score
istream & cath::index::filter::operator>>(istream              &arg_istream,             ///< The istream from which to extract the filter_vs_full_score
                                          filter_vs_full_score &arg_filter_vs_full_score ///< The filter_vs_full_score to populate
                                          ) {
	// Grab the string
	string input_string;
	getline( arg_istream,  input_string );

	// Trim off any prefix/suffix at the start/end
	if ( starts_with( input_string, filter_vs_full_score::STRING_PREFIX ) && ends_with( input_string, filter_vs_full_score::STRING_SUFFIX ) ) {
		const size_t prefix_length        = filter_vs_full_score::STRING_PREFIX.length();
		const size_t prefix_suffix_length = filter_vs_full_score::STRING_SUFFIX.length() + prefix_length;
		input_string = input_string.substr(
			filter_vs_full_score::STRING_PREFIX.length(),
			input_string.length() - min( input_string.length(), prefix_suffix_length )
		);
	}

	// Trim off any whitespace at the start/end
	trim( input_string );

	// Split on space or comma characters
	const str_vec input_strings = split_build<str_vec>(
		input_string,
		is_any_of( ", " ),
		token_compress_on
	);

	// Check that there are two parts
	if ( input_strings.size() != 2 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot parse filter_vs_full_score from string \""
			+ input_string
			+ "\" that contains "
			+ to_string( input_strings.size() )
			+ " part(s), rather than 2"
		));
	}

	// Convert the two strings to the filter and full scores and use them
	// to construct a filter_vs_full_score
	const double filter_score = stod( input_strings[ 0 ] );
	const double full_score   = stod( input_strings[ 1 ] );
	arg_filter_vs_full_score  = filter_vs_full_score( filter_score, full_score );

	// If all has gone well, then return the istream
	return arg_istream;
}

/// \brief Simple insertion operator for filter_vs_full_score
///
/// \relates filter_vs_full_score
ostream & cath::index::filter::operator<<(ostream                    &arg_os,                  ///< The ostream to which the filter_vs_full_score should be output
                                          const filter_vs_full_score &arg_filter_vs_full_score ///< The filter_vs_full_score to output
                                          ) {
	ostringstream temp_ss;
	temp_ss << filter_vs_full_score::STRING_PREFIX;
	temp_ss << right << setw( 7 ) << arg_filter_vs_full_score.get_filter_score();
	temp_ss << ",";
	temp_ss << right << setw( 7 ) << arg_filter_vs_full_score.get_full_score();
	temp_ss << filter_vs_full_score::STRING_SUFFIX;
	arg_os << temp_ss.str();
	return arg_os;
}
