/// \file
/// \brief The string_options_block class definitions

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

#include "string_options_block.hpp"

#include "cath/common/clone/make_uptr_clone.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::boost::program_options::options_description;
using ::boost::program_options::variables_map;
using ::std::nullopt;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone() method to act as a virtual copy-ctor
///
/// This is a concrete definition of a virtual method that's pure in options_block
unique_ptr<options_block> string_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Provide a name for the options block, as used in the options description text
///
/// This is a concrete definition of a virtual method that's pure in options_block
string string_options_block::do_get_block_name() const {
	return the_string;
}

/// \brief Add the block's non-hidden options to the provided options_description
///
/// This string_options_block  nothing
///
/// This is a concrete definition of a virtual method that's pure in options_block
void string_options_block::do_add_visible_options_to_description(options_description &/*prm_desc*/,       ///< The options_description to which the options are added
                                                                 const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                 ) {
}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// \returns A string describing the conflict in the options or an empty string if there's none
str_opt string_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                ) const {
	return nullopt;
}

/// \brief Return all options names for this block
///
/// This block has no options so this just returns empty
str_vec string_options_block::do_get_all_options_names() const {
	return {};
}

/// \brief Ctor from the string that this options block should put in the options usage text
string_options_block::string_options_block(string prm_string ///< The string that this options block should put in the options usage text
                                           ) : the_string{ std::move( prm_string ) } {
}

/// \brief Getter for the string that this options block should put in the options usage text
const string & string_options_block::get_string() const {
	return the_string;
}
