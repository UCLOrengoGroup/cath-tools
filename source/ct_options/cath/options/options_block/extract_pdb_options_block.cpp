/// \file
/// \brief The extract_pdb_options_block class definitions

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

#include "extract_pdb_options_block.hpp"

#include <filesystem>

#include "cath/common/clone/make_uptr_clone.hpp"
// #include "exception/invalid_argument_exception.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std::literals::string_literals;

using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::filesystem::path;
using ::std::nullopt;
using ::std::string;
using ::std::unique_ptr;

const string extract_pdb_options_block::PO_INPUT_PDB_FILE  ( "input-pdb-file"  );
const string extract_pdb_options_block::PO_OUTPUT_PDB_FILE ( "output-pdb-file" );
const string extract_pdb_options_block::PO_REGIONS         ( "regions"         );

/// \brief A standard do_clone method
///
/// This is a concrete definition of a virtual method that's pure in options_block
unique_ptr<options_block> extract_pdb_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
///
/// This is a concrete definition of a virtual method that's pure in options_block
string extract_pdb_options_block::do_get_block_name() const {
	return "Extract";
}

/// \brief Add this block's options to the provided options_description
///
/// This is a concrete definition of a virtual method that's pure in options_block
void extract_pdb_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                      const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                      ) {
	const auto out_pdb_notifier = [&] (const path   &x) { output_pdb_file = x; };
	const auto regions_notifier = [&] (const domain &x) { regions         = x; };

	const string file_varname    = "<file>";
	const string regions_varname = "<regions>";

	prm_desc.add_options()
		(
			PO_OUTPUT_PDB_FILE.c_str(),
			value<path>()
				->notifier     ( out_pdb_notifier )
				->value_name   ( file_varname ),
			( "Permit success for a file " + file_varname  + " that has no ATOM records" ).c_str()
		)
		(
			( PO_REGIONS ).c_str(),
			value<domain>()
				->notifier     ( regions_notifier )
				->value_name   ( regions_varname ),
			( "Extract region(s) " + regions_varname + " from the structure.\n"
				+ "Format is: 251-348:B,408-416A:B"
				// + "(Put " + regions_varname + R"( in quotes to prevent the square brackets confusing your shell ("No match")))"
				).c_str()
		);
}

/// \brief Add this block's hidden options to the provided options_description
void extract_pdb_options_block::do_add_hidden_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                     const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                     ) {
	const string file_varname = "<file>";

	prm_desc.add_options()
		(
			PO_INPUT_PDB_FILE.c_str(),
			value<path>( &input_pdb_file)
				// ->notifier  ( in_pdb_notifier )
				->value_name( file_varname    ),
			"PDB file to read from"
		);
}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// This is a concrete definition of a virtual method that's pure in options_block
///
/// At present, this always accepts all options
str_opt extract_pdb_options_block::do_invalid_string(const variables_map &prm_variables_map ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                     ) const {
	// If there is no PDB file to extract then grumble
	//
	// (Best done here rather than via boost::program_options::typed_value::required() because
	//  that leads to an error message that's unclear for users that don't know about positional options
	//  being implemented via hidden options)
	if ( ! specifies_option( prm_variables_map, PO_INPUT_PDB_FILE ) ) {
		return "Must specify an input PDB file to extract."s;
	}

	// Otherwise return all OK
	return nullopt;
}

/// \brief Return all options names for this block
str_vec extract_pdb_options_block::do_get_all_options_names() const {
	return {
		extract_pdb_options_block::PO_INPUT_PDB_FILE,
		extract_pdb_options_block::PO_OUTPUT_PDB_FILE,
		extract_pdb_options_block::PO_REGIONS,
	};
}

// /// \brief Getter for the PDB file
// path extract_pdb_options_block::get_pdb_file() const {
// 	return pdb_file;
// }

// /// \brief Getter for whether to permit no atoms
// bool extract_pdb_options_block::get_permit_no_atoms() const {
// 	return permit_no_atoms;
// }

/// \brief TODOCUMENT
const path & extract_pdb_options_block::get_input_pdb_file() const {
	return input_pdb_file;
}

/// \brief TODOCUMENT
const path_opt & extract_pdb_options_block::get_output_pdb_file() const {
	return output_pdb_file;
}

/// \brief TODOCUMENT
const domain_opt & extract_pdb_options_block::get_regions() const {
	return regions;
}

