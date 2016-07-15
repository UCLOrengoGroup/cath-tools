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

#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>

#include "common/boost_addenda/range/adaptor/lexical_casted.h"
#include "common/type_aliases.h"
#include "exception/program_exception_wrapper.h"
#include "resolve_hits/hit_arch.h"
#include "resolve_hits/hit_list.h"
#include "resolve_hits/hit_resolver.h"
#include "resolve_hits/read_and_resolve_mgr.h"
#include "resolve_hits/res_arrow.h"
#include "resolve_hits/scored_hit_arch.h"

#include <chrono>

using namespace cath::common;
using namespace cath::rslv;

using boost::algorithm::join;
using boost::filesystem::exists;
using boost::filesystem::path;
using std::chrono::high_resolution_clock;
using std::cin;
using std::cout;
using std::string;

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
				std::cerr << "Usage: cath-resolve-hits hits_file_to_process\n (use 'cath-resolve-hits -'' to read from stdin)\n";
				return;
			}

			// If dash, then read input from stdin
			if ( argv[ 1 ] == string{ "-" } ) {
				// Read the hits from stdin and print the resolved hits to cout
				read_and_resolve_mgr the_read_and_resolve_mgr{ cout };
				read_hit_list_from_istream( the_read_and_resolve_mgr, cin );
			}

			// Grab the filename and check the file exists
			const path the_file{ argv[ 1 ] };
			if ( ! exists( the_file ) ) {
				std::cerr << "Error: no such file \"" << the_file << "\" exists\n";
				return;
			}

			// Read the hits and print the resolved hits to cout
			read_and_resolve_mgr the_read_and_resolve_mgr{ cout };
			read_hit_list_from_file( the_read_and_resolve_mgr, the_file );
		}
	};
}

/// \brief A main function for cath_resolve_hits that just calls run_program() on a cath_resolve_hits_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath::cath_resolve_hits_program_exception_wrapper().run_program( argc, argv );
}
