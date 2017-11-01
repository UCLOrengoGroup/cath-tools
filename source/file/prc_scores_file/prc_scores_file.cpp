/// \file
/// \brief The prc_scores_file class definitions

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

#include "prc_scores_file.hpp"

#include <boost/algorithm/string/trim.hpp>
#include <boost/log/trivial.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/remove_copy_if.hpp>

#include "common/algorithm/contains.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/file/open_fstream.hpp"
#include "common/hash/pair_hash.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "file/prc_scores_file/detail/prc_scores_line_parser.hpp"
#include "file/prc_scores_file/prc_scores_entry.hpp"

#include <fstream>
#include <iostream> // ***** TEMPORARY *****
#include <string>
#include <unordered_map>

using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::adaptors::filtered;
using boost::algorithm::trim_left;
using boost::filesystem::path;

/// \brief TODOCUMENT
prc_scores_entry_vec prc_scores_file::remove_duplicates(const prc_scores_entry_vec &arg_prc_scores_entries ///< TODOCUMENT
                                                        ) {
	prc_scores_entry_vec  results;

	unordered_map<str_str_pair, size_t, pair_hash> index_of_previously_seen;

	return transform_build<prc_scores_entry_vec>(
		indices( arg_prc_scores_entries.size() )
			| filtered(
				[&] (const size_t &x) {
					const auto &entry     = arg_prc_scores_entries[ x ];
					const auto &id1       = entry.get_name_1();
					const auto &id2       = entry.get_name_2();
					const auto &evalue    = entry.get_evalue();
					const auto  name_pair = make_pair( id1, id2 );

					const auto find_itr = index_of_previously_seen.find( name_pair );

					if ( find_itr == common::cend( index_of_previously_seen ) ) {
						index_of_previously_seen.emplace( name_pair, x );
						return true;
					}

					const auto prev_entry = arg_prc_scores_entries[ find_itr->second ];
//					if ( entry.get_hit_num() <= prev_entry.get_hit_num() ) {
//						BOOST_THROW_EXCEPTION(invalid_argument_exception(
//							"When parsing PRC results, found hit between " + id1 + " and " + id2 + " with hit number that isn't higher than for previous result"
//						));
//						BOOST_LOG_TRIVIAL( warning ) << "When parsing PRC results, found hit between " << id1 << " and " << id2 << " with hit number that isn't higher than for previous result";
//					}
					if ( evalue          <  prev_entry.get_evalue()  ) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"When parsing PRC results, found hit between " + id1 + " and " + id2 + " with evalue better than for previous result"
						));
					}
					if ( entry.get_simple()  >  prev_entry.get_simple()  ) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"When parsing PRC results, found hit between " + id1 + " and " + id2 + " with simple score better than for previous result "
							+ " " + std::to_string(          entry.get_simple() )
							+ " " + std::to_string( prev_entry.get_simple() )
						));
					}
					return false;
				}
			),
		[&] (const size_t &x) { return arg_prc_scores_entries[ x ]; }
	);
	return results;
}

/// \brief Parse a vector of prc_scores_entry objects from the specified istream
prc_scores_entry_vec prc_scores_file::parse_prc_scores_file(istream &arg_prc_scores_is ///< The istream of prc scores data from which to parse the prc_scores_entry objects
                                                            ) {
	string                 line_string;
	prc_scores_entry_vec   results;
	prc_scores_line_parser parser;
	while ( getline( arg_prc_scores_is, line_string ) ) {

		// If this line is neither empty nor a comment line (a comment line is a line with a '#' character as the first non-whitespace character)
		trim_left( line_string );
		if ( ! line_string.empty() && line_string.front() != '#' ) {
			results.push_back( parser.parse_line( line_string ) );
		}
	}
	return results;
}

/// \brief Parse a vector of prc_scores_entry objects from the specified string
prc_scores_entry_vec prc_scores_file::parse_prc_scores_file(const string &arg_prc_scores_str ///< The string of prc scores data from which to parse the prc_scores_entry objects
                                                            ) {
	istringstream input_ss{ arg_prc_scores_str };
	return parse_prc_scores_file( input_ss );
}

/// \brief TODOCUMENT
prc_scores_entry_vec prc_scores_file::parse_prc_scores_file(const path &arg_prc_scores_file ///< TODOCUMENT
                                                            ) {
	ifstream prc_scores_ifstream;
	open_ifstream( prc_scores_ifstream, arg_prc_scores_file );
	const auto prc_scores_data = prc_scores_file::parse_prc_scores_file( prc_scores_ifstream );
	prc_scores_ifstream.close();
	return prc_scores_data;
}


/// \brief Parse a vector of prc_scores_entry objects from the specified istream
prc_scores_entry_vec prc_scores_file::parse_prc_scores_file_fancy(istream &arg_prc_scores_is ///< The istream of prc scores data from which to parse the prc_scores_entry objects
                                                                  ) {
	return remove_duplicates( parse_prc_scores_file( arg_prc_scores_is ) );
}

/// \brief Parse a vector of prc_scores_entry objects from the specified string
prc_scores_entry_vec prc_scores_file::parse_prc_scores_file_fancy(const string &arg_prc_scores_str ///< The string of prc scores data from which to parse the prc_scores_entry objects
                                                                  ) {
	istringstream input_ss{ arg_prc_scores_str };
	return parse_prc_scores_file_fancy( input_ss );
}

/// \brief TODOCUMENT
prc_scores_entry_vec prc_scores_file::parse_prc_scores_file_fancy(const path &arg_prc_scores_file ///< TODOCUMENT
                                                                  ) {
	return remove_duplicates( parse_prc_scores_file( arg_prc_scores_file ) );
}
