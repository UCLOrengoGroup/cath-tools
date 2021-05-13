/// \file
/// \brief The pdb_input_options_block class definitions

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

#include "pdb_input_options_block.hpp"

#include <filesystem>

#include "cath/common/clone/make_uptr_clone.hpp"

using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::cath;

using ::boost::program_options::bool_switch;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::filesystem::path;
using ::std::nullopt;
using ::std::string;
using ::std::unique_ptr;

/// \brief The option name for the a list of PDB files that should be read
const string pdb_input_options_block::PO_PDB_INFILE     ( "pdb-infile"      );

/// \brief The option name for whether to read PDBs from stdin
const string pdb_input_options_block::PO_PDBS_FROM_STDIN( "pdbs-from-stdin" );

/// \brief A standard do_clone method
unique_ptr<options_block> pdb_input_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string pdb_input_options_block::do_get_block_name() const {
	return "PDB files source";
}

/// \brief Add this block's options to the provided options_description
void pdb_input_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                    const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                    ) {
	const string pdb_file_varname = "<pdbfile>";

	const auto pdb_infile_notifier = [&] (const path_vec &x) { the_pdb_input_spec.set_input_files    ( x ); };
	const auto stdin_read_notifier = [&] (const bool     &x) { the_pdb_input_spec.set_read_from_stdin( x ); };

	prm_desc.add_options()
		(
			PO_PDB_INFILE.c_str(),
			value<path_vec>()
				->value_name( pdb_file_varname        )
				->notifier  ( pdb_infile_notifier ),
			( "Read PDB from file " + pdb_file_varname + " (may be specified multiple times)" ).c_str()
		)
		(
			PO_PDBS_FROM_STDIN.c_str(),
			bool_switch()
				->notifier     ( stdin_read_notifier                     )
				->default_value( pdb_input_spec::DEFAULT_READ_FROM_STDIN ),
			"Read PDBs from stdin (separated by line: \"END   \")"
		);

	static_assert( ! pdb_input_spec::DEFAULT_READ_FROM_STDIN,
		"If pdb_input_spec::DEFAULT_READ_FROM_STDIN isn't false, it might mess up the bool switch in here" );
}

str_opt pdb_input_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                   ) const {
	for (const path &input_file : the_pdb_input_spec.get_input_files() ) {
		if ( ! is_acceptable_input_file( input_file ) ) {
			return "PDB input file " + input_file.string() + " is not a valid input file";
		}
	}
	return nullopt;
}

/// \brief Return all options names for this block
str_vec pdb_input_options_block::do_get_all_options_names() const {
	return {
		pdb_input_options_block::PO_PDB_INFILE,
		pdb_input_options_block::PO_PDBS_FROM_STDIN,
	};
}

/// \brief TODOCUMENT
const pdb_input_spec & pdb_input_options_block::get_pdb_input_spec() const {
	return the_pdb_input_spec;
}

/// \brief Get the number of pdb_acquirer objects implied by the pdb_input_spec in the specified pdb_input_options_block
///
/// \relates pdb_input_options_block
size_t cath::opts::get_num_acquirers(const pdb_input_options_block &prm_pdb_input_options_block ///< The pdb_input_options_block to query
                                     ) {
	return get_num_acquirers( prm_pdb_input_options_block.get_pdb_input_spec() );
}

