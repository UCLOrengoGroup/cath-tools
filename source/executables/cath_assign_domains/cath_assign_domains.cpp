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
#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>

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
//using namespace cath::opts;
using namespace std;

using boost::algorithm::is_space;
using boost::algorithm::token_compress_on;
using boost::filesystem::path;
using boost::is_space;

namespace cath {
	/// \brief TODOCUMENT
	vector<pair<path, path> > parse_ssap_and_prc_files_data(const path &arg_filename ///< TODOCUMENT
	                                                        ) {
		ifstream data_data_ifstream;
		open_ifstream( data_data_ifstream, arg_filename );

		vector<pair<path, path> > data;
		string line_string;
		while ( getline( data_data_ifstream, line_string ) ) {
			const auto line_parts = split_build<str_vec>( line_string, is_space(), token_compress_on );
			data.emplace_back( line_parts[ 0 ], line_parts[ 1 ] );
//			const auto line_parts[ 0 ], line_parts[ 1 ];
		}
		data_data_ifstream.close();

		return data;
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
			const path data_data_file = "/export/people/ucbctnl/ticket_914_data/data_data.txt";
			const path sf_of_dom_file = "/export/people/ucbctnl/ticket_914_data/superfamily_of_domain.txt";

			BOOST_LOG_TRIVIAL( warning ) << "About to parse_superfamily_of_domain";
			const auto sf_of_dom = parse_superfamily_of_domain  ( sf_of_dom_file );
			BOOST_LOG_TRIVIAL( warning ) << "About to parse_ssap_and_prc_files_data";
			const auto bob       = parse_ssap_and_prc_files_data( data_data_file );

//			boost::optional<std::reference_wrapper<const ssap_and_prc>> best_magic_function(const ssaps_and_prcs_of_query &);

			BOOST_LOG_TRIVIAL( warning ) << "About to loop";
			for (const auto &x : bob) {
				const path &ssap_file = x.first;
				const path &prc_file  = x.second;
				const auto the_ssaps_and_prcs = make_ssaps_and_prcs_of_query(
					ssap_scores_file::parse_ssap_scores_file_simple( ssap_file ),
					prc_scores_file::parse_prc_scores_file_fancy   ( prc_file  )
				);
				const auto best_result = first_result_if(
					the_ssaps_and_prcs,
					[] (const ssap_and_prc &x, const ssap_and_prc &y) {
						// Reverse inequality to put the highest magic_function values to the start
						return x.get_magic_function_score() > y.get_magic_function_score();
					},
					[] (const ssap_and_prc &) { return true; }
				);
				if ( best_result ) {
					cerr << "best result is : " << best_result->get() << "\n";
				}
				else {
					cerr << "No best result" << endl;
				}


			}

			cerr << "Hello world\n";
		}
	};
}

/// \brief A main function for cath_assign_domains that just calls run_program() on a cath_assign_domains_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath::cath_assign_domains_program_exception_wrapper().run_program( argc, argv );
}


