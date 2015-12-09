/// \file
/// \brief The cath_assign_domains_program_exception_wrapper definitions

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

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>

#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/file/open_fstream.h"
#include "common/type_aliases.h"
#include "exception/program_exception_wrapper.h"
#include "file/prc_scores_file/prc_scores_file.h"
#include "file/ssap_scores_file/ssap_scores_file.h"
#include "score/homcheck_tools/first_result_if.h"
#include "score/homcheck_tools/ssaps_and_prcs_of_query.h"
#include "score/homcheck_tools/superfamily_of_domain.h"

#include <fstream>

using namespace cath::common;
using namespace cath::file;
using namespace cath::homcheck;
using namespace std;

using boost::algorithm::contains;
using boost::algorithm::is_space;
using boost::algorithm::join;
using boost::algorithm::token_compress_on;
using boost::filesystem::path;
using boost::is_space;
using boost::lexical_cast;
using boost::log::trivial::info;
using boost::log::trivial::severity;
using boost::make_optional;

namespace cath {

	/// \brief Parse a file with each line containing two, whitespace-separated entries: the SSAP and PRC files
	vector<pair<path, path> > parse_ssap_and_prc_files_data(const path &arg_filename ///< The file to parse
	                                                        ) {
		ifstream data_data_ifstream;
		open_ifstream( data_data_ifstream, arg_filename );

		vector<pair<path, path> > data;
		string line_string;
		while ( getline( data_data_ifstream, line_string ) ) {
			const auto line_parts = split_build<str_vec>( line_string, is_space(), token_compress_on );
			data.emplace_back( line_parts[ 0 ], line_parts[ 1 ] );
		}
		data_data_ifstream.close();

		return data;
	}

	/// \brief The query ID extracted from either SSAP or PRC results, whichever is non-empty
	///        or boost::none if both are empty
	opt_str query_id_of_either(const ssap_scores_entry_vec &arg_ssaps, ///< The possibly-empty SSAP results
	                           const prc_scores_entry_vec  &arg_prcs   ///< The possibly-empty PRC  results
	                           ) {
		return ( ! arg_prcs.empty()                   ) ? make_optional( arg_prcs.front().get_name_1()  ) :
		       ( ! arg_ssaps.empty()                  ) ? make_optional( arg_ssaps.front().get_name_1() ) :
		                                                  boost::none;
	}

	/// \brief The query length extracted from either SSAP or PRC results, whichever is non-empty
	///        or boost::none if both are empty
	opt_size query_legnth_of_either(const ssap_scores_entry_vec &arg_ssaps, ///< The possibly-empty SSAP results
	                                const prc_scores_entry_vec  &arg_prcs   ///< The possibly-empty PRC  results
	                                ) {
		return ( ! arg_prcs.empty()                   ) ? make_optional( arg_prcs.front().get_length_1()  ) :
		       ( ! arg_ssaps.empty()                  ) ? make_optional( arg_ssaps.front().get_length_1() ) :
		                                                  boost::none;
	}

	/// \brief Get the best specified number of PRC hits to domains in CATH
	prc_scores_entry_vec best_n_prc_hits(const prc_scores_entry_vec  &arg_prcs,      ///< The PRC results for the query in question (the query must be entry 1 in all PRC results)
	                                     const superfamily_of_domain &arg_sf_of_dom, ///< The superfamily_of_domain for finding which matches are assigned
	                                     const size_t                &arg_n          ///< The maximum number of results to return
	                                     ) {
		return first_n_results_if(
			arg_prcs,
			[] (const prc_scores_entry &x, const prc_scores_entry &y) {
				// Keep less-than inequality to put the smallest PRC evalues to the start
				return x.get_evalue() < y.get_evalue();
			},
			[&] (const prc_scores_entry &x) {
				return arg_sf_of_dom.has_superfamily_of_domain( x.get_name_2() );
			},
			arg_n
		);
	}

	/// \brief Get the best specified number of SSAP hits to domains in CATH
	ssap_scores_entry_vec best_n_ssap_hits(const ssap_scores_entry_vec &arg_ssaps,     ///< The SSAP results for the query in question (the query must be entry 1 in all SSAP results)
	                                       const superfamily_of_domain &arg_sf_of_dom, ///< The superfamily_of_domain for finding which matches are assigned
	                                       const size_t                &arg_n          ///< The maximum number of results to return
	                                       ) {
		return first_n_results_if(
			arg_ssaps,
			[] (const ssap_scores_entry &x, const ssap_scores_entry &y) {
				// Reverse inequality to put the highest SSAP scores to the start
				return x.get_ssap_score() > y.get_ssap_score();
			},
			[&] (const ssap_scores_entry &x) {
				return arg_sf_of_dom.has_superfamily_of_domain( x.get_name_2() );
			},
			arg_n
		);
	}

	/// \brief Generate Trac Wiki string describing a new fold result
	///
	/// ||=  Query Domain  =||=  Query Length  =||=  Best SSAP Hits  =||=  Best PRC Hits  =||
	string string_for_new_fold_strings(const path                  &arg_ssap_file, ///< The file from which the SSAP results have been parsed
	                                   const path                  &arg_prc_file,  ///< The file from which the PRC  results have been parsed
	                                   const ssap_scores_entry_vec &arg_ssaps,     ///< The SSAP results
	                                   const prc_scores_entry_vec  &arg_prcs,      ///< The PRC  results
	                                   const superfamily_of_domain &arg_sf_of_dom  ///< The superfamily_of_domain information
	                                   ) {
		constexpr size_t NUM_BEST_HITS = 6;
		const auto query_id     = query_id_of_either    ( arg_ssaps, arg_prcs );
		const auto query_length = query_legnth_of_either( arg_ssaps, arg_prcs );
		if ( ! query_id || ! query_length ) {
			return "|| COMPLETELY EMPTY RESULTS FOR " + arg_ssap_file.string() + " and " + arg_prc_file.string() + "||";
		}

		// Generate string describing this query's best SSAP hits
		const auto best_ssap_strings = transform_build<str_vec>(
			best_n_ssap_hits( arg_ssaps, arg_sf_of_dom, NUM_BEST_HITS ),
			[&] (const ssap_scores_entry &x) {
				return x.get_name_2()
					+ "; SSAP:**" + ( boost::format( "%g") % x.get_ssap_score() ).str()
					+ "**; O/L:**"  + ( boost::format( "%g") % x.get_overlap_pc() ).str()
					+ "**; RMSD:**" + ( boost::format( "%g") % x.get_rmsd()       ).str()
					+ "**; S/F:**"  + arg_sf_of_dom.get_superfamily_of_domain( x.get_name_2() )
					+ "**";
			}
		);

		// Generate string describing this query's best PRC hits
		const auto best_prc_strings = transform_build<str_vec>(
			best_n_prc_hits( arg_prcs, arg_sf_of_dom, NUM_BEST_HITS ),
			[&] (const prc_scores_entry &x) {
				return x.get_name_2()
					+ "; PRC:**"  + ( boost::format( "%g") % x.get_evalue() ).str()
					+ "**; S/F:**"  + arg_sf_of_dom.get_superfamily_of_domain( x.get_name_2() )
					+ "**";
			}
		);

		// Generate a Trac Wiki string
		return "||"
			+ *query_id
			+ " [[BR]] homcheck:"
			+ *query_id
			+ " [[BR]] domchop:"
			+ query_id->substr( 0, 5 )
			+ " [[BR]] rasmol-struc:"
			+ *query_id
			+ "  ||  "
			+ to_string( *query_length )
			+ "  ||"
			+ join( best_ssap_strings, "[[BR]]" )
			+ "  ||"
			+ join( best_prc_strings, "[[BR]]" )
			+ "  ||";
	}

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to cath_align_refiner::refine()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class cath_assign_domains_program_exception_wrapper final : public program_exception_wrapper {
		virtual string do_get_program_name() const override final {
			return "cath-assign-domains";
		}

		/// \brief Parse the options and then pass them to cath_assign_domainsr::superpose()
		virtual void do_run_program(int /*argc*/, char * /*argv*/[]) override final {
			boost::log::core::get()->set_filter(
				severity >= info
			);

			const path data_data_file = "/data1/people/ucbctnl/ticket_914_data/data_data.txt";
			const path sf_of_dom_file = "/data1/people/ucbctnl/ticket_914_data/superfamily_of_domain.txt";

			BOOST_LOG_TRIVIAL( info ) << "About to parse superfamily_of_domain data";
			auto sf_of_dom = parse_superfamily_of_domain  ( sf_of_dom_file );

			BOOST_LOG_TRIVIAL( info ) << "About to parse ssap_and_prc_files data";
			const auto ssap_and_prc_files = parse_ssap_and_prc_files_data( data_data_file );

			str_vec new_fold_strings;
			for (const auto &ssap_and_prc_file_pair : ssap_and_prc_files) {
				const path &ssap_file          = ssap_and_prc_file_pair.first;
				const path &prc_file           = ssap_and_prc_file_pair.second;
				const auto the_ssaps           = ssap_scores_file::parse_ssap_scores_file_simple( ssap_file );
				const auto the_prcs            = prc_scores_file::parse_prc_scores_file_fancy   ( prc_file  );
				const auto the_ssaps_and_prcs  = make_ssaps_and_prcs_of_query  ( the_ssaps, the_prcs );
				const auto best_result_ref_opt = best_magic_function_assignable( the_ssaps_and_prcs, sf_of_dom );
				if ( best_result_ref_opt ) {
					const auto &best_result    = best_result_ref_opt->get();
					const string &query_id     = best_result.get_query_id();
					const string &match_id     = best_result.get_match_id();
					const string &match_sf     = sf_of_dom.get_superfamily_of_domain( match_id );
					const string  match_fold   = homcheck::detail::fold_of_superfamily_id( match_sf );
					const bool    created_sf   = contains( match_sf, superfamily_of_domain::NEW_SF_CORE_STRING );
					const string  action       = created_sf ? "ASSIGN_TO_RECENT_SUPERFAMILY"
					                                        : "ASSIGN_TO_EXISTG_SUPERFAMILY";
					const string  match_node   = created_sf ? match_fold
					                                      : match_sf;
					const double &ssap_score   = get_ssap_score     ( best_result );
					const double &ssap_overlap = get_ssap_overlap_pc( best_result );
					const double &ssap_rmsd    = get_ssap_rmsd      ( best_result );
					const double &prc_evalue   = get_prc_evalue     ( best_result );
					const double &magic_fn     = best_result.get_magic_function_score();

					cout << query_id     << " "
					     << action       << " "
					     << match_id     << " "
					     << match_node   << " \"Assigning "
					     << query_id     << " to "
					     << match_node   << " based on strong match with "
					     << match_id     << " (SSAP score: "
					     << ssap_score   << "; SSAP overlap: "
					     << ssap_overlap << "; SSAP RMSD: "
					     << ssap_rmsd    << "; PRC evalue: "
					     << prc_evalue   << "; CMF: "
					     << magic_fn     << ")\"\n";
				}
				else {
					const auto fold_match_ref_opt = best_fold_level_match( the_ssaps, sf_of_dom );
					if ( fold_match_ref_opt ) {
						const auto &fold_match     = fold_match_ref_opt->get();
						const string &query_id     = fold_match.get_name_1();
						const string &match_id     = fold_match.get_name_2();
						const string &match_sf     = sf_of_dom.get_superfamily_of_domain( match_id );
						const string  match_node   = homcheck::detail::fold_of_superfamily_id( match_sf );
						const string  action       = "CREATE_SF_IN_FOLD_AND_ASSIGN";
						const double &ssap_score   = fold_match.get_ssap_score();
						const double &ssap_overlap = fold_match.get_overlap_pc();
						const double &ssap_rmsd    = fold_match.get_rmsd();

						cout << query_id     << " "
						     << action       << " "
						     << match_id     << " "
						     << match_node   << " \"Assigning "
						     << query_id     << " to new superfamily in "
						     << match_node   << " based on moderate SSAP match with "
						     << match_id     << " (SSAP score: "
						     << ssap_score   << "; SSAP overlap: "
						     << ssap_overlap << "; SSAP RMSD: "
						     << ssap_rmsd    << ")\"\n";

						sf_of_dom.add_domain_in_new_sf_in_fold_of_domain( query_id, match_id );
					}
				}
			}
		}
	};
}

/// \brief A main function for cath_assign_domains that just calls run_program() on a cath_assign_domains_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath::cath_assign_domains_program_exception_wrapper().run_program( argc, argv );
}


