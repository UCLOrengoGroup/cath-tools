/// \file
/// \brief The hmmer_scores_file class definitions

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

#include "hmmer_scores_file.hpp"

#include <filesystem>
#include <fstream>
#include <iostream> // ***** TEMPORARY *****
#include <string>

#include <boost/algorithm/string/trim.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/remove_copy_if.hpp>

#include <spdlog/spdlog.h>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/file/hmmer_scores_file/hmmer_scores_entry.hpp"


using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::std;

using ::boost::adaptors::filtered;
using ::boost::algorithm::trim_left;
using ::std::filesystem::path;

/// \brief TODOCUMENT
hmmer_scores_entry_vec hmmer_scores_file::remove_duplicates(const hmmer_scores_entry_vec &prm_hmmer_scores_entries ///< TODOCUMENT
                                                            ) {
	hmmer_scores_entry_vec  results;

	str_str_pair_size_map index_of_previously_seen;

	return transform_build<hmmer_scores_entry_vec>(
		indices( prm_hmmer_scores_entries.size() )
			| filtered(
				[&] (const size_t &x) {
					const auto &entry     = prm_hmmer_scores_entries[ x ];
					const auto &id1       = entry.get_name_1();
					const auto &id2       = entry.get_name_2();
					const auto &evalue    = entry.get_full_sequence_evalue();
					const auto  name_pair = make_pair( id1, id2 );

					if ( ! contains( index_of_previously_seen, name_pair ) ) {
						index_of_previously_seen.emplace( name_pair, x );
						return true;
					}

					cerr << "Ooo - crazy things happening with HMMER parsing" << "\n";

					const auto prev_entry = prm_hmmer_scores_entries[ index_of_previously_seen.at( name_pair ) ];
//					if ( entry.get_hit_num() <= prev_entry.get_hit_num() ) {
//						BOOST_THROW_EXCEPTION(invalid_argument_exception(
//							"When parsing PRC results, found hit between " + id1 + " and " + id2 + " with hit number that isn't higher than for previous result"
//						));
//						::spdlog::warn(
//							"When parsing PRC results, found hit between {} and {} with hit number that isn't higher than for previous result",
//							id1,
//							id2
//						);
//					}
					if ( evalue          <  prev_entry.get_full_sequence_evalue()  ) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"When parsing HMMER results, found hit between " + id1 + " and " + id2 + " with evalue better than for previous result"
						));
					}
					if ( entry.get_full_sequence_score()  >  prev_entry.get_full_sequence_score()  ) {
						BOOST_THROW_EXCEPTION(invalid_argument_exception(
							"When parsing HMMER results, found hit between " + id1 + " and " + id2 + " with simple score better than for previous result "
							+ " " + std::to_string( entry.get_full_sequence_score() )
							+ " " + std::to_string( prev_entry.get_full_sequence_score() )
						));
					}
					return false;
				}
			),
		[&] (const size_t &x) { return prm_hmmer_scores_entries[ x ]; }
	);
	return results;
}

/// \brief Parse a vector of hmmer_scores_entry objects from the specified istream
hmmer_scores_entry_vec hmmer_scores_file::parse_hmmer_scores_file(istream                   &prm_hmmer_scores_is,    ///< The istream of prc scores data from which to parse the hmmer_scores_entry objects
                                                                  const hmmer_name_handling &prm_hmmer_name_handling ///< TODOCUMENT
                                                                  ) {
	string               line_string;
	hmmer_scores_entry_vec results;
	while ( getline( prm_hmmer_scores_is, line_string ) ) {

		// If this line is neither empty nor a comment line (a comment line is a line with a '#' character as the first non-whitespace character)
		trim_left( line_string );
		if ( ! line_string.empty() && line_string.front() != '#' ) {
			results.push_back( hmmer_scores_entry_from_line( line_string, prm_hmmer_name_handling ) );
		}
	}
	return results;
}

/// \brief TODOCUMENT
hmmer_scores_entry_vec hmmer_scores_file::parse_hmmer_scores_file(const path                &prm_hmmer_scores_file,  ///< TODOCUMENT
                                                                  const hmmer_name_handling &prm_hmmer_name_handling ///< TODOCUMENT
                                                                  ) {
	ifstream hmmer_scores_ifstream = open_ifstream( prm_hmmer_scores_file );
	const auto hmmer_scores_data = hmmer_scores_file::parse_hmmer_scores_file( hmmer_scores_ifstream, prm_hmmer_name_handling );
	hmmer_scores_ifstream.close();
	return hmmer_scores_data;
}


///// \brief Parse a vector of hmmer_scores_entry objects from the specified istream
//hmmer_scores_entry_vec hmmer_scores_file::parse_hmmer_scores_file_fancy(istream &prm_hmmer_scores_is ///< The istream of prc scores data from which to parse the hmmer_scores_entry objects
//                                                                        ) {
//	return remove_duplicates( parse_hmmer_scores_file( prm_hmmer_scores_is ) );
//}

///// \brief TODOCUMENT
//hmmer_scores_entry_vec hmmer_scores_file::parse_hmmer_scores_file_fancy(const path &prm_hmmer_scores_file ///< TODOCUMENT
//                                                                        ) {
//	return remove_duplicates( parse_hmmer_scores_file( prm_hmmer_scores_file ) );
//}
