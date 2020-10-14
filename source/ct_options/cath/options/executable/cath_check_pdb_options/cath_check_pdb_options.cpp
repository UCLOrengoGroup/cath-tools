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

#include "cath_check_pdb_options.hpp"

#include <boost/program_options.hpp>
#include <boost/shared_array.hpp>

#include "cath/common/argc_argv_faker.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

#include <fstream>

using namespace cath;
using namespace cath::common;
using namespace cath::opts;

using boost::filesystem::path;
using boost::none;
using boost::program_options::positional_options_description;
using std::ifstream;
using std::string;

/// \brief The name of the program that uses this executable_options
const string cath_check_pdb_options::PROGRAM_NAME("check-pdb");

/// \brief Get the name of the program that uses this executable_options
string cath_check_pdb_options::do_get_program_name() const {
	return PROGRAM_NAME;
}

/// \brief TODOCUMENT
positional_options_description cath_check_pdb_options::get_positional_options() {
	positional_options_description positionals;
	positionals.add( check_pdb_options_block::PO_PDB_FILE.c_str(), 1 );
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
str_opt cath_check_pdb_options::do_get_error_or_help_string() const {
	return none;
}

/// \brief Get a string to prepend to the standard help
string cath_check_pdb_options::do_get_help_prefix_string() const {
	return "Usage: " + PROGRAM_NAME + " pdb_file\n\n"
		+ get_overview_string();
}

/// \brief Get a string to append to the standard help (just empty here)
string cath_check_pdb_options::do_get_help_suffix_string() const {
	return "";
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string cath_check_pdb_options::do_get_overview_string() const {
	return "Check a PDB file for some potential problems";
}

/// \brief Ctor for cath_check_pdb_options
cath_check_pdb_options::cath_check_pdb_options() {
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


