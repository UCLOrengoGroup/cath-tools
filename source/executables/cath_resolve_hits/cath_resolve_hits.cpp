/// \file
/// \brief The cath_resolve_hits_program_exception_wrapper definitions

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

//#include <iostream> // ***** TEMPORARY ****
//using std::cerr;    // ***** TEMPORARY ****

#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>

#include "common/boost_addenda/range/adaptor/lexical_casted.h"
#include "common/type_aliases.h"
#include "exception/program_exception_wrapper.h"
#include "resolve_hits/hit_arch.h"
#include "resolve_hits/hit_list.h"
#include "resolve_hits/hit_resolver.h"
#include "resolve_hits/res_arrow.h"
#include "resolve_hits/scored_hit_arch.h"

#include <chrono>

using namespace cath::common;
using namespace cath::rslv;

using boost::algorithm::join;
using boost::filesystem::exists;
using boost::filesystem::path;
using std::string;
using std::chrono::high_resolution_clock;

namespace cath {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to cath_hits_resolver::resolve()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class cath_resolve_hits_program_exception_wrapper final : public program_exception_wrapper {
		virtual string do_get_program_name() const override final {
			return "cath-resolve-hits";
		}

		/// \brief Parse the options and then pass them to cath_resolve_hits::superpose()
		virtual void do_run_program(int argc, char * argv[]) override final {
			if ( argc != 2 ) {
				std::cerr << "Usage: cath-resolve-hits hits_file_to_process\n";
				return;
			}

			const path the_file{ argv[ 1 ] };
			if ( ! exists( the_file ) ) {
				std::cerr << "Error: no such file \"" << the_file << "\" exists\n";
				return;
			}

			const auto      start_time  = high_resolution_clock::now();
			const auto      the_hits    = read_hit_list_from_file( the_file );
			const auto      loaded_time = high_resolution_clock::now();
			const auto      best_result = resolve_hits( the_hits );
			const auto      finish_time = high_resolution_clock::now();
			const resscr_t &best_score  = best_result.get_score();
			std::cerr << "Input file                   : " << the_file                 << "\n";
			std::cerr << "Number of input hits         : " << the_hits.size()          << "\n";
			std::cerr << "Maximum stop residue number  : " << get_max_stop( the_hits ) << "\n";
			std::cerr << "Time spent [parsing hits   ] : " << durn_to_seconds_string( loaded_time - start_time  ) << "\n";
			std::cerr << "Time spent [finding optimum] : " << durn_to_seconds_string( finish_time - loaded_time ) << "\n";
			std::cerr << "Time spent [total          ] : " << durn_to_seconds_string( finish_time - start_time ) << "\n";
			std::cerr << "Final score                  : " << best_score << std::endl;
			std::cout << to_string( best_result.get_arch(), hit_output_format::JON );

			// // std::cerr << "file is : " << file << "\n";
			// // std::cerr << "Exists  : " << std::boolalpha << exists( file ) << "\n";
			// // return;
			// const path_vec the_files = {
			// };

			// for (const auto &the_file : the_files) {
			// 	const auto      start_time  = high_resolution_clock::now();
			// 	const auto      the_hits    = read_hit_list_from_file( the_file );
			// 	const auto      best_result = resolve_hits( the_hits );
			// 	const resscr_t &best_score  = best_result.get_score();
			// 	// const hit_arch &best_arch   = best_result.get_arch();
			// 	std::cerr
			// 		<< "Durn: "
			// 		<< durn_to_seconds_string( high_resolution_clock::now() - start_time )
			// 		<< "; size: "
			// 		<< the_hits.size()
			// 		<< "; seq_length: "
			// 		<< get_max_stop( the_hits )
			// 		<< "; result_score: "
			// 		<< best_score
			// 		<< "; file: "
			// 		<< the_file
			// 		// << "\n"
			// 		// << best_arch
			// 		<< "\n";
			// }

//			for (const auto &hit_a : the_hits) {
//				if ( hit_a.is_discontig() ) {
//					for (const auto &hit_b : the_hits) {
//						if ( hit_a.get_label() != hit_b.get_label() && hit_b.is_discontig() ) {
//							if ( any_interaction( hit_a, hit_b ) && ! cath::rslv::hits_overlap( hit_a, hit_b ) ) {
//								std::cerr << "any_interaction: "
//									<< ( second_right_intersperses_first( hit_a, hit_b ) ? "***" : "   " )
//									<< "\t" << get_segments_string( hit_a ) << "\t" << get_segments_string( hit_b ) << "\n";
//									std::cerr << "\t" << join(
//										rslv::detail::hit_resolver{ the_hits }.get_arrows_before_starts_of_doms_right_interspersed_with_all_of( { hit_a } ) | lexical_casted<string>(),
//										"\n\t"
//									) << "\n";
////								}
//							}
//
//						}
//					}
//				}
//			}

			// BOOST_LOG_TRIVIAL( warning )
			// 	<< "final result is "
			// 	<< best_score
			// 	<< "\n"
			// 	<< best_arch
			// 	<< "\n";
		}
	};
}

/// \brief A main function for cath_resolve_hits that just calls run_program() on a cath_resolve_hits_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath::cath_resolve_hits_program_exception_wrapper().run_program( argc, argv );
}
