/// \file
/// \brief The pdb_input_options_block class definitions

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

#include "pdb_input_options_block.h"

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/optional.hpp>

#include "common/clone/make_uptr_clone.h"
#include "options/acquirer/pdbs_acquirer/istream_pdbs_acquirer.h"
#include "options/acquirer/pdbs_acquirer/file_list_pdbs_acquirer.h"

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::assign::ptr_push_back;
using boost::none;
using boost::ptr_vector;

const string pdb_input_options_block::PO_PDB_INFILE     ( "pdb-infile"      );
const string pdb_input_options_block::PO_PDBS_FROM_STDIN( "pdbs-from-stdin" );

/// \brief A standard do_clone method.
unique_ptr<options_block> pdb_input_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string pdb_input_options_block::do_get_block_name() const {
	return "PDB files source";
}

/// \brief Add this block's options to the provided options_description
void pdb_input_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                    ) {
	arg_desc.add_options()
		(PO_PDB_INFILE.c_str(),      value<path_vec >(&input_files),                      "Read PDB from file arg (specify multiple times)"      )
		(PO_PDBS_FROM_STDIN.c_str(), bool_switch(&read_from_stdin)->default_value(false), "Read PDBs from stdin (separated by line: \"END   \")" );
}

opt_str pdb_input_options_block::do_invalid_string() const {
	const path_vec input_files = get_input_files_cref();
	for (const path &input_file : input_files) {
		if (!is_acceptable_input_file(input_file)) {
			return "PDB input file " + input_file.string() + " is not a valid input file";
		}
	}
	return none;
}

/// \brief TODOCUMENT
const path_vec & pdb_input_options_block::get_input_files_cref() const {
	return input_files;
}

/// \brief TODOCUMENT
bool pdb_input_options_block::get_read_from_stdin() const {
	return read_from_stdin;
}

/// \brief TODOCUMENT
ptr_vector<pdbs_acquirer> pdb_input_options_block::get_pdbs_acquirers() const {
	ptr_vector<pdbs_acquirer> pdb_acquirers;
	if ( get_read_from_stdin() ) {
		ptr_push_back<istream_pdbs_acquirer>( pdb_acquirers )();
	}
	if ( ! get_input_files_cref().empty() ) {
		ptr_push_back<file_list_pdbs_acquirer>( pdb_acquirers )( get_input_files_cref() );
	}
	return pdb_acquirers;
}

