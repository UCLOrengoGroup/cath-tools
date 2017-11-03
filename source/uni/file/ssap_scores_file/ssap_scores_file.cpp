/// \file
/// \brief The ssap_scores_file class definitions

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

#include "ssap_scores_file.hpp"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "common/algorithm/contains.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/file/open_fstream.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "file/ssap_scores_file/ssap_scores_entry.hpp"

#include <fstream>
#include <map>

using namespace boost::algorithm;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::algorithm::is_any_of;
using boost::algorithm::token_compress_on;
using boost::filesystem::path;
using boost::lexical_cast;

/// \brief Parse a vector of ssap_scores_entry objects from the specified istream
ssap_scores_entry_vec ssap_scores_file::parse_ssap_scores_file_simple(istream &arg_ssap_scores_is ///< The istream of SSAP scores data from which to parse the ssap_scores_entry objects
                                                                      ) {
	string                line_string;
	ssap_scores_entry_vec results;
	while ( getline( arg_ssap_scores_is, line_string ) ) {

		// If this line is neither empty nor a comment line (a comment line is a line with a '#' character as the first non-whitespace character)
		trim_left( line_string );
		if ( ! line_string.empty() && line_string.front() != '#' ) {
			results.push_back( ssap_scores_entry_from_line( line_string ) );
		}
	}
	return results;
}

/// \brief Parse a vector of ssap_scores_entry objects from the specified string
ssap_scores_entry_vec ssap_scores_file::parse_ssap_scores_file_simple(const string &arg_ssap_scores_string ///< TODOCUMENT
                                                                      ) {
	istringstream input_ss{ arg_ssap_scores_string };
	return parse_ssap_scores_file_simple( input_ss );
}

/// \brief Parse a vector of ssap_scores_entry objects from the specified file
ssap_scores_entry_vec ssap_scores_file::parse_ssap_scores_file_simple(const path &arg_ssap_scores_file ///< The SSAP scores file from which to parse the ssap_scores_entry objects
                                                                      ) {
	ifstream ssap_scores_ifstream;
	open_ifstream( ssap_scores_ifstream, arg_ssap_scores_file );
	const auto ssap_scores_data = ssap_scores_file::parse_ssap_scores_file_simple( ssap_scores_ifstream );
	ssap_scores_ifstream.close();
	return ssap_scores_data;
}

/// \brief TODOCUMENT
pair<str_vec, size_size_doub_tpl_vec> ssap_scores_file::parse_ssap_scores_file(istream &arg_ssap_scores_is ///< TODOCUMENT
                                                                               ) {
	const auto entries = parse_ssap_scores_file_simple( arg_ssap_scores_is );

	str_vec               id1s;
	str_vec               id2s;
	str_str_pair_doub_map score_of_ids;
	for (const auto &entry : entries) {
		// Grab the two IDs and the SSAP score and store them
		const auto &id1 = entry.get_name_1();
		const auto &id2 = entry.get_name_2();
		id1s.push_back( id1 );
		id2s.push_back( id2 );
		score_of_ids[ make_pair( id1, id2 ) ] = entry.get_ssap_score();
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
	size_size_doub_tpl_vec score_of_index_pair;
	for (const auto &score_data : score_of_ids) {
		const string &id1        = score_data.first.first;
		const string &id2        = score_data.first.second;
		const double &ssap_score = score_data.second;
		const size_t  index_1    = index_of_id[ id1 ];
		const size_t  index_2    = index_of_id[ id2 ];
		score_of_index_pair.emplace_back( index_1, index_2, ssap_score );

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
pair<str_vec, size_size_doub_tpl_vec> ssap_scores_file::parse_ssap_scores_file(const path &arg_ssap_scores_file ///< TODOCUMENT
                                                                               ) {
	ifstream ssap_scores_ifstream;
	open_ifstream( ssap_scores_ifstream, arg_ssap_scores_file );
	const pair<str_vec, size_size_doub_tpl_vec> ssap_scores_data = ssap_scores_file::parse_ssap_scores_file(ssap_scores_ifstream);
	ssap_scores_ifstream.close();
	return ssap_scores_data;
}

/// \brief Create arbitrary positive/negative values for all pair in the specified ssap_scores_entries for testing purposes
///
/// This alternates negative and positive
str_str_pair_bool_map cath::file::make_arbitrary_is_positive_data(const ssap_scores_entry_vec &arg_ssap_scores_entries ///< The arg_ssap_scores_entries from which to extract pairs to be made negative or positive
                                                                  ) {
	return transform_build<str_str_pair_bool_map>(
		indices( arg_ssap_scores_entries.size() ),
		[&] (const size_t &x) {
			const auto &the_ssap_scores_entry = arg_ssap_scores_entries[ x ];
			return make_pair(
				make_pair( the_ssap_scores_entry.get_name_1(), the_ssap_scores_entry.get_name_2() ),
				( x % 2 == 0 )
			);
		}
	);
}
