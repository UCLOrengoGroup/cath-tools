/// \file
/// \brief The hmmer_scores_file class definitions

#include "hmmer_scores_file.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/log/trivial.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/remove_copy_if.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/contains.h"
#include "common/algorithm/transform_build.h"
#include "common/file/open_fstream.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "file/hmmer_scores_file/hmmer_scores_entry.h"

#include <fstream>
#include <iostream> // ***** TEMPORARY *****
#include <string>

using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::adaptors::filtered;
using boost::algorithm::trim_left;
using boost::filesystem::path;
using boost::irange;

/// \brief TODOCUMENT
hmmer_scores_entry_vec hmmer_scores_file::remove_duplicates(const hmmer_scores_entry_vec &arg_hmmer_scores_entries ///< TODOCUMENT
                                                            ) {
	hmmer_scores_entry_vec  results;

	str_str_pair_size_map index_of_previously_seen;

	return transform_build<hmmer_scores_entry_vec>(
		irange( 0_z, arg_hmmer_scores_entries.size() )
			| filtered(
				[&] (const size_t &x) {
					const auto &entry     = arg_hmmer_scores_entries[ x ];
					const auto &id1       = entry.get_name_1();
					const auto &id2       = entry.get_name_2();
					const auto &evalue    = entry.get_full_sequence_evalue();
					const auto  name_pair = make_pair( id1, id2 );

					if ( ! contains( index_of_previously_seen, name_pair ) ) {
						index_of_previously_seen.emplace( name_pair, x );
						return true;
					}

					cerr << "Ooo - crazy things happening with HMMER parsing" << "\n";

					const auto prev_entry = arg_hmmer_scores_entries[ index_of_previously_seen.at( name_pair ) ];
//					if ( entry.get_hit_num() <= prev_entry.get_hit_num() ) {
//						BOOST_THROW_EXCEPTION(invalid_argument_exception(
//							"When parsing PRC results, found hit between " + id1 + " and " + id2 + " with hit number that isn't higher than for previous result"
//						));
//						BOOST_LOG_TRIVIAL( warning ) << "When parsing PRC results, found hit between " << id1 << " and " << id2 << " with hit number that isn't higher than for previous result";
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
		[&] (const size_t &x) { return arg_hmmer_scores_entries[ x ]; }
	);
	return results;
}

/// \brief Parse a vector of hmmer_scores_entry objects from the specified istream
hmmer_scores_entry_vec hmmer_scores_file::parse_hmmer_scores_file(istream                   &arg_hmmer_scores_is,    ///< The istream of prc scores data from which to parse the hmmer_scores_entry objects
                                                                  const hmmer_name_handling &arg_hmmer_name_handling ///< TODOCUMENT
                                                                  ) {
	string               line_string;
	hmmer_scores_entry_vec results;
	while ( getline( arg_hmmer_scores_is, line_string ) ) {

		// If this line is neither empty nor a comment line (a comment line is a line with a '#' character as the first non-whitespace character)
		trim_left( line_string );
		if ( ! line_string.empty() && line_string.front() != '#' ) {
			results.push_back( hmmer_scores_entry_from_line( line_string, arg_hmmer_name_handling ) );
		}
	}
	return results;
}

/// \brief TODOCUMENT
hmmer_scores_entry_vec hmmer_scores_file::parse_hmmer_scores_file(const path                &arg_hmmer_scores_file,  ///< TODOCUMENT
                                                                  const hmmer_name_handling &arg_hmmer_name_handling ///< TODOCUMENT
                                                                  ) {
	ifstream hmmer_scores_ifstream;
	open_ifstream( hmmer_scores_ifstream, arg_hmmer_scores_file );
	const auto hmmer_scores_data = hmmer_scores_file::parse_hmmer_scores_file( hmmer_scores_ifstream, arg_hmmer_name_handling );
	hmmer_scores_ifstream.close();
	return hmmer_scores_data;
}


///// \brief Parse a vector of hmmer_scores_entry objects from the specified istream
//hmmer_scores_entry_vec hmmer_scores_file::parse_hmmer_scores_file_fancy(istream &arg_hmmer_scores_is ///< The istream of prc scores data from which to parse the hmmer_scores_entry objects
//                                                                        ) {
//	return remove_duplicates( parse_hmmer_scores_file( arg_hmmer_scores_is ) );
//}

///// \brief TODOCUMENT
//hmmer_scores_entry_vec hmmer_scores_file::parse_hmmer_scores_file_fancy(const path &arg_hmmer_scores_file ///< TODOCUMENT
//                                                                        ) {
//	return remove_duplicates( parse_hmmer_scores_file( arg_hmmer_scores_file ) );
//}
