/// \file
/// \brief The label_pair_is_positive class definitions

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

#include "cath/score/score_classification/label_pair_is_positive/label_pair_is_positive.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>

#include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/file/open_fstream.hpp"

#include <fstream>

using namespace ::cath::common;
using namespace ::cath::score;
using namespace ::std;

using ::boost::algorithm::is_space;
using ::boost::algorithm::trim_left;
using ::boost::filesystem::path;
using ::boost::lexical_cast;
using ::boost::token_compress_on;

/// \brief Ctor from a map of pairs of labels to is_positive bools
label_pair_is_positive::label_pair_is_positive(str_str_pair_bool_map prm_positive_of_label_pairs ///< A map of pairs of labels to is_positive bools from which to construct this label_pair_is_positive
                                               ) : all_pairs { std::move( prm_positive_of_label_pairs ) } {

}

/// \brief Query whether a given pair is positive
///
/// \pre An entry should be stored for the specified pair else an invalid_argument_exception will be thrown
bool label_pair_is_positive::is_positive(const string &prm_label_a, ///< The first  label to query
                                         const string &prm_label_b  ///< The second label to query
                                         ) const {
	const auto find_itr = all_pairs.find( make_pair( prm_label_a, prm_label_b ) );
	if ( find_itr == common::cend( all_pairs ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to find label pair \""
			+ prm_label_a
			+ "\", \""
			+ prm_label_b
			+ "\" in label_pair_is_positive"
		));
	}
	return find_itr->second;
}

/// \brief Parse a label_pair_is_positive from an istream that has one line per entry, each with two labels and a true/false or 1/0
///
/// \relates label_pair_is_positive
label_pair_is_positive cath::score::make_label_pair_is_positive(istream &prm_input_text ///< The istream of data from which to parse the label_pair_is_positive
                                                                ) {
	string                line_string;
	str_str_pair_bool_map results;
	while ( getline( prm_input_text, line_string ) ) {
		// If this line isn't empty...
		trim_left( line_string );
		if ( ! line_string.empty() ) {
			const auto line_parts = split_build<str_vec>( line_string, is_space(), token_compress_on );
			if ( line_parts.size() != 3 ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to parse label_pair_is_positive from input with line that doesn't contain three entries"));
			}
			const bool positive = lexical_cast<bool>( line_parts[ 2 ] );
			results.insert( make_pair( make_pair( line_parts[ 0 ], line_parts[ 1 ] ), positive ) );
		}
	}
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return label_pair_is_positive{ results };
}

/// \brief Parse a label_pair_is_positive from a file that has one line per entry, each with two labels and a true/false or 1/0
///
/// \relates label_pair_is_positive
label_pair_is_positive cath::score::make_label_pair_is_positive(const path &prm_input_file ///< The file of data from which to parse the label_pair_is_positive
                                                                ) {
	ifstream ssap_scores_ifstream;
	open_ifstream( ssap_scores_ifstream, prm_input_file );
	const auto the_label_pair_is_positive = make_label_pair_is_positive(ssap_scores_ifstream);
	ssap_scores_ifstream.close();
	return the_label_pair_is_positive;
}
