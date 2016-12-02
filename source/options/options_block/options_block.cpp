/// \file
/// \brief The options_block class definitions

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

#include "options_block.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional.hpp>

#include "common/clone/check_uptr_clone_against_this.hpp"

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::filesystem::is_empty;
using boost::none;
using boost::numeric_cast;


/// \brief A string to use to separate (valid values and their descriptions) from each other
const string options_block::SUB_DESC_SEPARATOR = "\n   ";

/// \brief The string to use to separate valid values' names from their descriptions
const string options_block::SUB_DESC_PAIR_SEPARATOR = " - ";

/// \brief TODOCUMENT
void options_block::do_add_hidden_options_to_description(options_description &/*arg_desc*/ ///< TODOCUMENT
                                                         ) {
}

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<options_block> options_block::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief A method that uses the concrete class's methods to construct an options description
///        (using the specified line_length) of the hidden options
options_description options_block::get_all_options_description(const size_t &arg_line_length ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                               ) {
	options_description desc{ do_get_block_name(), numeric_cast<unsigned int>( arg_line_length ) };
	do_add_visible_options_to_description( desc );
	do_add_hidden_options_to_description ( desc );
	return desc;
}

/// \brief A method that uses the concrete class's methods to construct an options description
///        (using the specified line_length) of the visible options
///
/// The Boost program_options documentation seems to be a bit out of sync with the actual code.
/// In particular, line_length isn't (well) documented.
options_description options_block::get_visible_options_description(const size_t &arg_line_length ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                   ) {
	const string block_name = do_get_block_name();
	options_description desc( block_name, numeric_cast<unsigned int>( arg_line_length ) );
	do_add_visible_options_to_description( desc );
	return desc;
}

/// \brief A method that uses the concrete class's methods to construct an options description
///        (using the specified line_length) of the hidden options
///
/// The Boost program_options documentation seems to be a bit out of sync with the actual code.
/// In particular, line_length isn't (well) documented.
options_description options_block::get_hidden_options_description() {
	const string block_name = do_get_block_name();
	options_description desc( block_name + " [hidden options]" );
	do_add_hidden_options_to_description( desc );
	return desc;
}

/// \brief An NVI pass-through to a method that returns a string describing any problems or "" if none
str_opt options_block::invalid_string(const variables_map &arg_variables_map ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                      ) const {
	return do_invalid_string( arg_variables_map );
}


/// \brief TODOCUMENT
bool options_block::is_acceptable_output_file(const path &arg_output_file ///< TODOCUMENT
                                              ) {
	if (exists(arg_output_file) && !is_regular_file(arg_output_file)) {
		return false;
	}
	if (arg_output_file.has_parent_path() && !is_directory(arg_output_file.parent_path())) {
		return false;
	}
	return true;
}

/// \brief TODOCUMENT
bool options_block::is_acceptable_input_file(const path &arg_output_file, ///< TODOCUMENT
                                             const bool &arg_allow_empty  ///< TODOCUMENT
                                             ) {
	using boost::filesystem::is_empty;

	if (!exists(arg_output_file)) {
		return false;
	}
	if (!arg_allow_empty && is_empty(arg_output_file)) {
		return false;
	}
	return true;
}

/// \brief TODOCUMENT
bool options_block::is_acceptable_output_dir(const path &arg_output_dir ///< TODOCUMENT
                                             ) {
	if (!exists(arg_output_dir) && !is_directory(arg_output_dir)) {
		return false;
	}
	return true;
}

/// \brief TODOCUMENT
bool options_block::is_acceptable_input_dir(const path &arg_input_dir ///< TODOCUMENT
                                            ) {
	if (!exists(arg_input_dir) && !is_directory(arg_input_dir)) {
		return false;
	}
	return true;
}

/// \brief TODOCUMENT
bool options_block::is_acceptable_executable(const path &arg_output_file ///< TODOCUMENT
                                             ) {
	if (arg_output_file.filename() == ".." || arg_output_file.filename() == ".") {
//		cerr << "Here0" << endl;
		return false;
	}
	if (arg_output_file.has_root_directory() || arg_output_file.has_root_name()) {
//		cerr << "Here1" << endl;
		if (!exists(arg_output_file) && !is_regular_file(arg_output_file)) {
//			cerr << "Here2" << endl;
			return false;
		}
	}
	return true;
}
