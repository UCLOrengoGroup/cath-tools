/// \file
/// \brief The cath_check_pdb_options class definitions

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

#include "cath_check_pdb_options.h"

#include <boost/program_options.hpp>
#include <boost/shared_array.hpp>

#include "common/argc_argv_faker.h"
#include "common/file/open_fstream.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"

#include <fstream>

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace std;

const string cath_check_pdb_options::PROGRAM_NAME("check-pdb");

/// TODOCUMENT
string cath_check_pdb_options::do_get_program_name() const {
	return PROGRAM_NAME;
}

/// TODOCUMENT
positional_options_description cath_check_pdb_options::get_positional_options() {
	positional_options_description positionals;
	positionals.add(check_pdb_options_block::PO_PDB_FILE.c_str(), -1);
	return positionals;
}

/// \brief Review all specified options and return a string containing any errors or a help string
///        (possibly using a description of all visible options)
///
/// This is a concrete definition of a virtual method that's pure in executable_options
///
/// This should only be called by executable_options, as the last step of the parse_options()
/// method, after all real parsing has completed.
///
/// \pre The options should have been parsed
///
/// \returns Any error/help string arising from the newly specified options
///          or an empty string if there aren't any
string cath_check_pdb_options::do_update_error_or_help_string(const options_description &arg_visible_program_options ///< The full options_description of visible options
                                                              ) const {
	// If help was requested, then provide it
	if (the_misc_options_block.get_help()) {
		return the_misc_options_block.get_help_string(arg_visible_program_options, get_help_prefix_string(), "");
	}

	// If version information was requested, then provide it
	if (the_misc_options_block.get_version()) {
		return the_misc_options_block.get_version_string( get_program_name(), get_overview_string() );
	}

	// If there is no PDB file to check then grumble
	const path pdb_file(get_pdb_file());
	if (pdb_file.empty()) {
		return "Must specify a PDB file to check.";
	}

	// Check the PDB file is a valid input file
	if (!options_block::is_acceptable_input_file(pdb_file)) {
		return "No such valid, non-empty PDB file \"" + pdb_file.string() + "\".";
	}

	// Otherwise, no problems have been detected so return an empty string
	return "";
}

string cath_check_pdb_options::get_help_prefix_string() {
	return "Usage: " + PROGRAM_NAME + " pdb_file\n\n"
		+ get_overview_string();
}

/// \brief Ctor for cath_check_pdb_options
cath_check_pdb_options::cath_check_pdb_options() {
	super::add_options_block( the_misc_options_block      );
	super::add_options_block( the_check_pdb_options_block );
}

/// \brief Getter for the PDB file
path cath_check_pdb_options::get_pdb_file() const {
	return the_check_pdb_options_block.get_pdb_file();
}

/// \brief Getter for whether to permit no ATOMs
bool cath_check_pdb_options::get_permit_no_atoms() const {
	return the_check_pdb_options_block.get_permit_no_atoms();
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string cath_check_pdb_options::get_overview_string() {
	return "Check a PDB file for some potential problems";
}

/// \brief Check that the PDB file is OK and throw an invalid_argument_exception if not
///
/// \returns Nothing
void cath::opts::check_pdb_file(const path &arg_pdb_file,       ///< The PDB file to check
                                const bool &arg_permit_no_atoms ///< Whether to permit no ATOM records
                                ) {
	// Check the PDB file is a valid input file
	if (!options_block::is_acceptable_input_file(arg_pdb_file)) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("No such valid, non-empty PDB file \"" + arg_pdb_file.string() + "\"."));
	}

	// Open an ifstream on the PDB files
	ifstream pdb_istream;
	open_ifstream(pdb_istream, arg_pdb_file);

	// Attempt to read the PDB file (and let any exceptions propagate out)
	const pdb newly_read_pdb = read_pdb_file( pdb_istream );
	pdb_istream.close();

	// If there were no ATOM records and that isn't allowed, then throw an exception
	// (which will be caught just below)
	if (!arg_permit_no_atoms && newly_read_pdb.get_num_atoms() <= 0) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("PDB file \"" + arg_pdb_file.string() + "\" did not contain any valid ATOM records"));
	}
}

