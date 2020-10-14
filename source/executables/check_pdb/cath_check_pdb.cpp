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

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/logger.hpp"
#include "cath/common/program_exception_wrapper.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/options/executable/cath_check_pdb_options/cath_check_pdb_options.hpp"

#include <fstream>

using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;

using ::boost::filesystem::path;
using ::std::cerr;
using ::std::cout;
using ::std::ifstream;
using ::std::string;

namespace cath {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to cath_superposer::superpose()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class cath_check_pdb_program_exception_wrapper final : public program_exception_wrapper {
		string do_get_program_name() const final {
			return "check-pdb";
		}

		/// \brief Parse the options and then pass them to cath_superposer::superpose()
		void do_run_program(int argc, char * argv[]) final {
			const auto the_cath_check_pdb_options = make_and_parse_options<cath_check_pdb_options>( argc, argv );

			// If the options are invalid or specify to do_nothing, then just output string return
			const auto &error_or_help_string = the_cath_check_pdb_options.get_error_or_help_string();
			if ( error_or_help_string ) {
				cout << *error_or_help_string;
				return;
			}

			const path pdb_file        = the_cath_check_pdb_options.get_pdb_file();
			const bool permit_no_atoms = the_cath_check_pdb_options.get_permit_no_atoms();
			try {
				check_pdb_file(pdb_file, permit_no_atoms);
			}
			// Catch any problems in reading the file
			catch (const invalid_argument_exception &prm_exception) {
				logger::log_and_exit(
					logger::return_code::MALFORMED_PDB_FILE,
					"Unable to parse PDB file \"" + pdb_file.string() + "\". Error was:\n" + prm_exception.what()
				);
			}

			// If no problems were caught, output that the PDB file was parsed successfully
			cerr << "PDB file " << pdb_file << " parsed successfully\n";
		}
	public:

		/// \brief Check that the PDB file is OK and throw an invalid_argument_exception if not
		///
		/// \returns Nothing
		static void check_pdb_file(const path &prm_pdb_file,       ///< The PDB file to check
		                           const bool &prm_permit_no_atoms ///< Whether to permit no ATOM records
		                           ) {
			// Check the PDB file is a valid input file
			if ( ! options_block::is_acceptable_input_file( prm_pdb_file ) ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("No such valid, non-empty PDB file \"" + prm_pdb_file.string() + "\"."));
			}

			// Open an ifstream on the PDB files
			ifstream pdb_istream;
			open_ifstream(pdb_istream, prm_pdb_file);

			// Attempt to read the PDB file (and let any exceptions propagate out)
			const pdb newly_read_pdb = read_pdb_file( pdb_istream );
			pdb_istream.close();

			// If there were no ATOM records and that isn't allowed, then throw an exception
			// (which will be caught just below)
			if (!prm_permit_no_atoms && newly_read_pdb.get_num_atoms() <= 0) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("PDB file \"" + prm_pdb_file.string() + "\" did not contain any valid ATOM records"));
			}
		}
	};
} // namespace cath

/// \brief A main function for cath_check_pdb that just calls run_program() on a cath_check_pdb_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath::cath_check_pdb_program_exception_wrapper().run_program( argc, argv );
}
