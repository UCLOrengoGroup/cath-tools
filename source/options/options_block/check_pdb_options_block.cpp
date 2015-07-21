/// \file
/// \brief The check_pdb_options_block class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "check_pdb_options_block.h"

#include <boost/optional.hpp>

#include "common/clone/make_uptr_clone.h"
#include "exception/invalid_argument_exception.h"

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::none;

const string check_pdb_options_block::PO_PDB_FILE ( "pdb-file"        );
const string check_pdb_options_block::PO_PERMIT   ( "permit-no-atoms" );

/// \brief A standard do_clone method
///
/// This is a concrete definition of a virtual method that's pure in options_block
unique_ptr<options_block> check_pdb_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
///
/// This is a concrete definition of a virtual method that's pure in options_block
string check_pdb_options_block::do_get_block_name() const {
	return "Check PDB";
}

/// \brief Add this block's options to the provided options_description
///
/// This is a concrete definition of a virtual method that's pure in options_block
void check_pdb_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                    ) {
	arg_desc.add_options()
		(PO_PDB_FILE.c_str(),  value<path>(&pdb_file),                              "PDB file to check"                                  )
		(PO_PERMIT.c_str(),    bool_switch(&permit_no_atoms)->default_value(false), "Permit success for a file that has no ATOM records" );
}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// This is a concrete definition of a virtual method that's pure in options_block
///
/// At present, this always accepts all options
opt_str check_pdb_options_block::do_invalid_string() const {
	return none;
}

/// \brief Construct a check_pdb_options_block from a map from option name to a pair of description and help message
check_pdb_options_block::check_pdb_options_block() : permit_no_atoms(false) {
}

/// \brief Getter for the PDB file
path check_pdb_options_block::get_pdb_file() const {
	return pdb_file;
}

/// \brief Getter for whether to permit no atoms
bool check_pdb_options_block::get_permit_no_atoms() const {
	return permit_no_atoms;
}
