/// \file
/// \brief The cath_check_pdb_program_exception_wrapper definitions

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

#include "common/logger.h"
#include "exception/invalid_argument_exception.h"
#include "exception/program_exception_wrapper.h"
#include "options/executable/cath_check_pdb_options/cath_check_pdb_options.h"

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

namespace cath {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to cath_superposer::superpose()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class cath_check_pdb_program_exception_wrapper final : public program_exception_wrapper {
		virtual string do_get_program_name() const override final {
			return "check-pdb";
		}

		/// \brief Parse the options and then pass them to cath_superposer::superpose()
		virtual void do_run_program(int argc, char * argv[]) override final {
			const auto the_cath_check_pdb_options = make_and_parse_options<cath_check_pdb_options>( argc, argv );

			// If the options are invalid or specify to do_nothing, then just output string return
			const string error_or_help_string = the_cath_check_pdb_options.get_error_or_help_string();
			if (!error_or_help_string.empty()) {
				cout << error_or_help_string << endl;
				return;
			}

			const path pdb_file        = the_cath_check_pdb_options.get_pdb_file();
			const bool permit_no_atoms = the_cath_check_pdb_options.get_permit_no_atoms();
			try {
				check_pdb_file(pdb_file, permit_no_atoms);
			}
			// Catch any problems in reading the file
			catch (const invalid_argument_exception &arg_exception) {
				logger::log_and_exit(
					logger::return_code::MALFORMED_PDB_FILE,
					"Unable to parse PDB file \"" + pdb_file.string() + "\". Error was:\n" + arg_exception.what()
				);
			}

			// If no problems were caught, output that the PDB file was parsed successfully
			cerr << "PDB file " << pdb_file << " parsed successfully" << endl;
		}
	};
}

/// \brief A main function for cath_check_pdb that just calls run_program() on a cath_check_pdb_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath::cath_check_pdb_program_exception_wrapper().run_program( argc, argv );
}
