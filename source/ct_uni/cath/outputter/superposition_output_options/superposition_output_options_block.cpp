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

#include "superposition_output_options_block.hpp"

#include <filesystem>
#include <iostream>

#include <fmt/core.h>

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/outputter/superposition_outputter/json_file_superposition_outputter.hpp"
#include "cath/outputter/superposition_outputter/ostream_superposition_outputter.hpp"
#include "cath/outputter/superposition_outputter/pdb_file_superposition_outputter.hpp"
#include "cath/outputter/superposition_outputter/pdb_files_superposition_outputter.hpp"
#include "cath/outputter/superposition_outputter/pymol_file_superposition_outputter.hpp"
#include "cath/outputter/superposition_outputter/pymol_view_superposition_outputter.hpp"
#include "cath/outputter/superposition_outputter/superposition_outputter_list.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::cath::sup;

using ::boost::program_options::bool_switch;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::filesystem::path;
using ::std::nullopt;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<options_block> superposition_output_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string superposition_output_options_block::do_get_block_name() const {
	return "Superposition output";
}

/// \brief Add this block's options to the provided options_description
void superposition_output_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                               const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                               ) {
	prm_desc.add_options()
		(string( PO_SUP_FILE ).c_str(),          value<path>(&sup_to_pdb_file),                                     "Write the superposed structures to a single PDB file arg, separated using faked chain codes"  )
		(string( PO_SUP_FILES_DIR ).c_str(),     value<path>(&sup_to_pdb_files_dir),                                "Write the superposed structures to separate PDB files in directory arg"                       )
		(string( PO_SUP_TO_STDOUT ).c_str(),     bool_switch(&sup_to_stdout)->default_value(false),                 "Print the superposed structures to stdout, separated using faked chain codes"                 )
		(string( PO_SUP_TO_PYMOL ).c_str(),      bool_switch(&sup_to_pymol )->default_value(false),                 "Start up PyMOL for viewing the superposition"                                                 )
		(string( PO_PYMOL_PROGRAM ).c_str(),     value<path>(&pymol_program)->default_value(DEFAULT_PYMOL_PROGRAM), "Use arg as the PyMOL executable for viewing; may optionally include the full path"            )
		(string( PO_SUP_TO_PYMOL_FILE ).c_str(), value<path>(&sup_to_pymol_file),                                   "Write the superposition to a PyMOL script arg\n(Recommended filename extension: .pml)"        )
		(string( PO_SUP_TO_JSON_FILE ).c_str(),  value<path>(&json_file),                                           "Write the superposition to JSON superposition file\n(Recommended filename extension: .sup_json)"        );
}

str_opt superposition_output_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                              ) const {
	if (!get_sup_to_pdb_file().empty() && !is_acceptable_output_file(get_sup_to_pdb_file())) {
		return "Not a valid superposition pdb output file :\"" + get_sup_to_pdb_file().string() + "\"";
	}
	if (!get_sup_to_pdb_files_dir().empty() && !is_acceptable_output_dir(get_sup_to_pdb_files_dir())) {
		return "Not a valid superposition pdb output dir :\"" + get_sup_to_pdb_files_dir().string() + "\"";
	}
	if (get_pymol_program() != DEFAULT_PYMOL_PROGRAM && !get_sup_to_pymol()) {
		return ::fmt::format( "Cannot specify PyMOL executable if not using --{}", PO_SUP_TO_PYMOL );
	}
	if (!is_acceptable_executable(get_pymol_program())) {
		return "Not a valid PyMOL executable :\"" + get_pymol_program().string() + "\"";
	}
	if (!get_sup_to_pymol_file().empty() && !is_acceptable_output_file(get_sup_to_pymol_file())) {
		return "Not a valid superposition PyMOL output file:\"" + get_sup_to_pymol_file().string() + "\"";
	}
	if (!get_json_file().empty() && !is_acceptable_output_file(get_json_file())) {
		return "Not a valid superposition JSON output file:\"" + get_json_file().string() + "\"";
	}

	return nullopt;
}

/// \brief Return all options names for this block
str_view_vec superposition_output_options_block::do_get_all_options_names() const {
	return {
		PO_SUP_FILE,
		PO_SUP_FILES_DIR,
		PO_SUP_TO_STDOUT,
		PO_SUP_TO_PYMOL,
		PO_PYMOL_PROGRAM,
		PO_SUP_TO_PYMOL_FILE,
		PO_SUP_TO_JSON_FILE,
	};
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

/// TODOCUMENT
path superposition_output_options_block::get_json_file() const {
	return json_file;
}

/// Get a list of superposition_outputter_list that will perform the outputs specified in this options_block
superposition_outputter_list superposition_output_options_block::get_superposition_outputters(const display_spec               &prm_display_spec,          ///< The specification of how to display (colour) the structures
                                                                                              const superposition_content_spec &prm_content_spec,          ///< The specification of what should be included in the superposition
                                                                                              const default_supn_outputter     &prm_default_supn_outputter ///< What superposition outputter (if any) should be provided if none is explicitly specified
                                                                                              ) const {
	superposition_outputter_list superposition_outputters;
	if ( ! get_sup_to_pdb_file().empty() ) {
		superposition_outputters.push_back( pdb_file_superposition_outputter   { get_sup_to_pdb_file(),                        prm_content_spec } );
	}
	if ( ! get_sup_to_pdb_files_dir().empty() ) {
		superposition_outputters.push_back( pdb_files_superposition_outputter  { get_sup_to_pdb_files_dir(),                   prm_content_spec } );
	}
	if ( get_sup_to_stdout() ) {
		superposition_outputters.push_back( ostream_superposition_outputter    {                                               prm_content_spec } );
	}
	if ( get_sup_to_pymol() ) {
		superposition_outputters.push_back( pymol_view_superposition_outputter { get_pymol_program(),        prm_display_spec, prm_content_spec } );
	}
	if ( ! get_sup_to_pymol_file().empty() ) {
		superposition_outputters.push_back( pymol_file_superposition_outputter { get_sup_to_pymol_file(),    prm_display_spec, prm_content_spec } );
	}
	if ( ! get_json_file().empty() ) {
		superposition_outputters.push_back( json_file_superposition_outputter  { get_json_file()                                                } );
	}

	// If there no superposition outputters have been specified,
	// add any defaults as requested
	if ( superposition_outputters.empty() ) {
		switch ( prm_default_supn_outputter ) {
			case ( default_supn_outputter::PYMOL ): {
				superposition_outputters.push_back(
				  pymol_view_superposition_outputter{ get_pymol_program(), prm_display_spec, prm_content_spec } );
				break;
			}
			case ( default_supn_outputter::NONE ): {
				break;
			}
		}
	}

	return superposition_outputters;
}

/// \brief TODOCUMENT
bool superposition_output_options_block::outputs_to_stdout() const {
	return get_sup_to_stdout();
}
