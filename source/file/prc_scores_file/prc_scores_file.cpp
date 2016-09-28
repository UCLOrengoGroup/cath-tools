/// \file
/// \brief The prc_scores_file class definitions

#include "prc_scores_file.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/log/trivial.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/remove_copy_if.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/contains.h"
#include "common/algorithm/transform_build.h"
#include "common/cpp14/cbegin_cend.h"
#include "common/file/open_fstream.h"
#include "common/hash/pair_hash.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "file/prc_scores_file/detail/prc_scores_line_parser.h"
#include "file/prc_scores_file/prc_scores_entry.h"

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
using boost::irange;

/// \brief TODOCUMENT
prc_scores_entry_vec prc_scores_file::remove_duplicates(const prc_scores_entry_vec &arg_prc_scores_entries ///< TODOCUMENT
                                                        ) {
	prc_scores_entry_vec  results;

	unordered_map<str_str_pair, size_t, pair_hash> index_of_previously_seen;

	return transform_build<prc_scores_entry_vec>(
		irange( 0_z, arg_prc_scores_entries.size() )
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
