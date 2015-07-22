/// \file
/// \brief The ssap_scores_file class definitions

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

#include "ssap_scores_file.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "common/algorithm/contains.h"
#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/file/open_fstream.h"
#include "exception/invalid_argument_exception.h"

#include <fstream>
#include <map>

using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::algorithm::is_any_of;
using boost::algorithm::token_compress_on;
using boost::lexical_cast;

/// \brief TODOCUMENT
pair<str_vec, size_size_pair_doub_map> ssap_scores_file::parse_ssap_scores_file(istream &arg_ssap_scores_is ///< TODOCUMENT
                                                                                ) {
	using str_str_pair = pair<string, string>;
	using str_str_pair_doub_map = map<str_str_pair, double>;

	string                line_string;
	str_vec               id1s;
	str_vec               id2s;
	str_str_pair_doub_map score_of_ids;
	while (getline(arg_ssap_scores_is, line_string)) {

		// If this line is neither empty nor a comment line (a comment line is a line with a '#' character as the first non-whitespace character)
		trim_left(line_string);
		if ( ! line_string.empty() && line_string.front() != '#' ) {

			// Split the line into parts
			const str_vec line_parts = split_build<str_vec>( line_string, is_any_of( " " ), token_compress_on );
			if ( line_parts.size() < 5 ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("SSAP scores line contains fewer than five entries but the SSAP score should be the fifth entry"));
			}

			// Grab the two ids and the SSAP score and store them
			const string id1        = line_parts[ 0 ];
			const string id2        = line_parts[ 1 ];
			const double ssap_score = stod( line_parts[ 4 ] );
			id1s.push_back( id1 );
			id2s.push_back( id2 );
			score_of_ids[ make_pair( id1, id2 ) ] = ssap_score;
		}
	}

	// Form a list of unique IDs ordered like the id1s and then the id2s
	str_vec ids;
	map<string, size_t> index_of_id;
	const str_vec_vec id1s_and_id2s = { id1s, id2s};
	for (const str_vec &id1s_or_id2s : id1s_and_id2s) {
		for (const string &id : id1s_or_id2s) {
			if ( ! contains( index_of_id, id ) ) {
//				pair<string, str_vec::size_type> bob = make_pair(id, ids.size());
				index_of_id[id] = ids.size();
				ids.push_back(id);
			}
		}
	}

	// Build a map of the scores based on the indices
	size_size_pair_doub_map score_of_index_pair;
	for (const str_str_pair_doub_map::value_type &score_data : score_of_ids) {
		const string &id1        = score_data.first.first;
		const string &id2        = score_data.first.second;
		const double &ssap_score = score_data.second;
		const size_t  index_1    = index_of_id[id1];
		const size_t  index_2    = index_of_id[id2];
		score_of_index_pair[make_pair(index_1, index_2)] = ssap_score;

		// TEMPORARY DEBUG STATEMENTS
//		cerr << "id1 is " << id1;
//		cerr << ", id2 is " << id2;
//		cerr << ", index_1 is " << index_1;
//		cerr << ", index_2 is " << index_2;
//		cerr << ", ssap_score is " << ssap_score;
//		cerr << endl;
	}

	// Return the result of all this hard work
	return make_pair(ids, score_of_index_pair);
}

/// \brief TODOCUMENT
pair<str_vec, size_size_pair_doub_map> ssap_scores_file::parse_ssap_scores_file(const path &arg_ssap_scores_file ///< TODOCUMENT
                                                                                ) {
	ifstream ssap_scores_ifstream;
	open_ifstream( ssap_scores_ifstream, arg_ssap_scores_file );
	const pair<str_vec, size_size_pair_doub_map> ssap_scores_data = ssap_scores_file::parse_ssap_scores_file(ssap_scores_ifstream);
	ssap_scores_ifstream.close();
	return ssap_scores_data;
}