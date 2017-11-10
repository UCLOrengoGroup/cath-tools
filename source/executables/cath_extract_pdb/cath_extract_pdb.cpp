/// \file
/// \brief The cath_extract_pdb definitions

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

#include <boost/filesystem.hpp>

#include "chopping/chopping_type_aliases.hpp"
#include "chopping/region/region.hpp"
#include "common/file/ofstream_list.hpp"
#include "common/program_exception_wrapper.hpp"
#include "file/pdb/pdb.hpp"
#include "options/executable/cath_extract_pdb_options/cath_extract_pdb_options.hpp"

using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;

using boost::filesystem::path;
using std::cout;
using std::string;

namespace cath {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to cath_superposer::superpose()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class cath_extract_pdb_program_exception_wrapper final : public program_exception_wrapper {
		string do_get_program_name() const final {
			return "cath-extract-pdb";
		}

		/// \brief Parse the options and then pass them to cath_superposer::superpose()
		void do_run_program(int    argc,
		                    char * argv[]
		                    ) final {
			const auto the_opts = make_and_parse_options<cath_extract_pdb_options>( argc, argv );

			// If the options are invalid or specify to do_nothing, then just output string return
			const auto &error_or_help_string = the_opts.get_error_or_help_string();
			if ( error_or_help_string ) {
				cout << *error_or_help_string;
				return;
			}

			const auto out_file_opt = get_output_pdb_file( the_opts );
			ofstream_list the_ofstreams{ std::cout };
			const path_vec paths =
				out_file_opt
				? path_vec{ { *out_file_opt } }
				: path_vec{ { path{ "-" }   } };
			auto ostreams = the_ofstreams.open_ofstreams( paths );

			write_pdb_file(
				ostreams[ 0 ],
				read_pdb_file( get_input_pdb_file( the_opts ) ).set_post_ter_residues( {} ),
				get_regions_opt( get_regions( the_opts ) )
			);
		}
	};
} // namespace cath

/// \brief A main function for cath_extract_pdb that just calls run_program() on a cath_extract_pdb_program_exception_wrapper
int main(int    argc,
         char * argv[] ) {
	return cath::cath_extract_pdb_program_exception_wrapper()
		.run_program( argc, argv );
}
