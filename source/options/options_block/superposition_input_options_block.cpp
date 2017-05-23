/// \file
/// \brief The superposition_input_options_block class definitions

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

#include "superposition_input_options_block.hpp"

#include "common/clone/make_uptr_clone.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::opts;

using boost::filesystem::path;
using boost::none;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief The option name for a file from which to read a JSON superposition
const string superposition_input_options_block::PO_JSON_SUP_INFILE{ "json-sup-infile" };

/// \brief A standard do_clone method
unique_ptr<options_block> superposition_input_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string superposition_input_options_block::do_get_block_name() const {
	return "Superposition source";
}

/// \brief Add this block's options to the provided options_description
void superposition_input_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                              ) {
	const string file_varname { "<file>" };

	const auto json_sup_infile_notifier = [&] (const path &x) { json_sup_infile = x; };

	arg_desc.add_options()
		(
			PO_JSON_SUP_INFILE.c_str(),
			value<path>()
				->value_name( file_varname             )
				->notifier  ( json_sup_infile_notifier ),
			( "Read superposition from file " + file_varname ).c_str()
		);
}

/// \brief Return a string describing any problems with the current configuration of the block
str_opt superposition_input_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                             ) const {
	return none;
}

/// \brief Return all options names for this block
str_vec superposition_input_options_block::do_get_all_options_names() const {
	return {
		superposition_input_options_block::PO_JSON_SUP_INFILE,
	};
}

/// \brief Getter for the optional file from which to read a JSON superposition
const path_opt & superposition_input_options_block::get_json_sup_infile() const {
	return json_sup_infile;
}
