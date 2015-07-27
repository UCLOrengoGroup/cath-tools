/// \file
/// \brief The superposition_output_options_block class definitions

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

#include "superposition_output_options_block.h"

#include <boost/optional.hpp>

#include "common/clone/make_uptr_clone.h"
#include "options/outputter/superposition_outputter/ostream_superposition_outputter.h"
#include "options/outputter/superposition_outputter/pdb_file_superposition_outputter.h"
#include "options/outputter/superposition_outputter/pdb_files_superposition_outputter.h"
#include "options/outputter/superposition_outputter/pymol_file_superposition_outputter.h"
#include "options/outputter/superposition_outputter/pymol_view_superposition_outputter.h"
#include "options/outputter/superposition_outputter/superposition_outputter_list.h"

#include <iostream>

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::none;

const string superposition_output_options_block::PO_SUP_FILE(          "sup-to-pdb-file"      );
const string superposition_output_options_block::PO_SUP_FILES_DIR(     "sup-to-pdb-files-dir" );
const string superposition_output_options_block::PO_SUP_TO_STDOUT(     "sup-to-stdout"        );
const string superposition_output_options_block::PO_SUP_TO_PYMOL(      "sup-to-pymol"         );
const string superposition_output_options_block::PO_PYMOL_PROGRAM(     "pymol-program"        );
const string superposition_output_options_block::PO_SUP_TO_PYMOL_FILE( "sup-to-pymol-file"    );

const string superposition_output_options_block::DEFAULT_PYMOL_PROGRAM("pymol");

/// \brief A standard do_clone method.
unique_ptr<options_block> superposition_output_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string superposition_output_options_block::do_get_block_name() const {
	return "Superposition output";
}

/// \brief Add this block's options to the provided options_description
void superposition_output_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                               ) {
	arg_desc.add_options()
		(PO_SUP_FILE.c_str(),          value<path>(&sup_to_pdb_file),                                     "Write the superposed structures to a single PDB file arg, separated using faked chain codes"  )
		(PO_SUP_FILES_DIR.c_str(),     value<path>(&sup_to_pdb_files_dir),                                "Write the superposed structures to separate PDB files in directory arg"                       )
		(PO_SUP_TO_STDOUT.c_str(),     bool_switch(&sup_to_stdout)->default_value(false),                 "Print the superposed structures to stdout, separated using faked chain codes"                 )
		(PO_SUP_TO_PYMOL.c_str(),      bool_switch(&sup_to_pymol )->default_value(false),                 "Start up PyMOL for viewing the superposition"                                                 )
		(PO_PYMOL_PROGRAM.c_str(),     value<path>(&pymol_program)->default_value(DEFAULT_PYMOL_PROGRAM), "Use arg as the PyMOL executable for viewing; may optionally include the full path"            )
		(PO_SUP_TO_PYMOL_FILE.c_str(), value<path>(&sup_to_pymol_file),                                   "Write the superposition to a PyMOL script arg\n(Recommended filename extension: .pml)"        );
}

opt_str superposition_output_options_block::do_invalid_string() const {
	if (!get_sup_to_pdb_file().empty() && !is_acceptable_output_file(get_sup_to_pdb_file())) {
		return "Not a valid superposition pdb output file :\"" + get_sup_to_pdb_file().string() + "\"";
	}
	if (!get_sup_to_pdb_files_dir().empty() && !is_acceptable_output_dir(get_sup_to_pdb_files_dir())) {
		return "Not a valid superposition pdb output dir :\"" + get_sup_to_pdb_files_dir().string() + "\"";
	}
	if (get_pymol_program() != DEFAULT_PYMOL_PROGRAM && !get_sup_to_pymol()) {
		return "Cannot specify PyMOL executable if not using --" + PO_SUP_TO_PYMOL + "";
	}
	if (!is_acceptable_executable(get_pymol_program())) {
		return "Not a valid PyMOL executable :\"" + get_pymol_program().string() + "\"";
	}
	if (!get_sup_to_pymol_file().empty() && !is_acceptable_output_file(get_sup_to_pymol_file())) {
		return "Not a valid superposition PyMOL output file:\"" + get_sup_to_pymol_file().string() + "\"";
	}

	return none;
}

/// TODOCUMENT
path superposition_output_options_block::get_sup_to_pdb_file() const {
	return sup_to_pdb_file;
}

/// TODOCUMENT
path superposition_output_options_block::get_sup_to_pdb_files_dir() const {
	return sup_to_pdb_files_dir;
}

/// TODOCUMENT
bool superposition_output_options_block::get_sup_to_stdout() const {
	return sup_to_stdout;
}

/// TODOCUMENT
bool superposition_output_options_block::get_sup_to_pymol() const {
	return sup_to_pymol;
}

/// TODOCUMENT
path superposition_output_options_block::get_pymol_program() const {
	return pymol_program;
}

/// TODOCUMENT
path superposition_output_options_block::get_sup_to_pymol_file() const {
	return sup_to_pymol_file;
}

superposition_outputter_list superposition_output_options_block::get_superposition_outputters(const display_spec &arg_display_spec ///< TODOCUMENT
                                                                                              ) const {
	superposition_outputter_list superposition_outputters;
	if ( ! get_sup_to_pdb_file().empty() ) {
		superposition_outputters.push_back( pdb_file_superposition_outputter ( get_sup_to_pdb_file() ) );
	}
	if ( ! get_sup_to_pdb_files_dir().empty() ) {
		superposition_outputters.push_back( pdb_files_superposition_outputter( get_sup_to_pdb_files_dir() ) );
	}
	if ( get_sup_to_stdout() ) {
		superposition_outputters.push_back( ostream_superposition_outputter() );
	}
	if ( get_sup_to_pymol() ) {
		superposition_outputters.push_back( pymol_view_superposition_outputter( get_pymol_program(), arg_display_spec ) );
	}
	if ( ! get_sup_to_pymol_file().empty() ) {
		superposition_outputters.push_back( pymol_file_superposition_outputter( get_sup_to_pymol_file(), arg_display_spec ) );
	}

	return superposition_outputters;
}

/// \brief TODOCUMENT
bool superposition_output_options_block::outputs_to_stdout() const {
	return get_sup_to_stdout();
}
