/// \file
/// \brief The detail_help_options_block class definitions

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

#include "detail_help_options_block.hpp"

#include <boost/optional.hpp>
#include <boost/range/adaptor/map.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::opts;

using boost::adaptors::map_keys;
using boost::none;
using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method
///
/// This is a concrete definition of a virtual method that's pure in options_block
unique_ptr<options_block> detail_help_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
///
/// This is a concrete definition of a virtual method that's pure in options_block
string detail_help_options_block::do_get_block_name() const {
	return "Detailed help";
}

/// \brief Add this block's options to the provided options_description
///
/// This is a concrete definition of a virtual method that's pure in options_block
void detail_help_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                      ) {
	// Foreach option requested
	for (const str_str_str_pair_pair &option_help_pair : desc_and_help_of_option_name) {
		// Grab the name, description and a reference to a (possibly newly created) bool value in the values map
		const string &option_name = option_help_pair.first;
		const string &description = option_help_pair.second.first;
		bool         &value_bool  = values[option_name];

		// Add this option to arg_desc
		arg_desc.add_options()
			(option_name.c_str(), bool_switch(&value_bool)->default_value(false), description.c_str());
	}
}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// This is a concrete definition of a virtual method that's pure in options_block
///
/// At present, this always accepts all options
str_opt detail_help_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                     ) const {
	return none;
}

/// \brief Return all options names for this block
str_vec detail_help_options_block::do_get_all_options_names() const {
	return copy_build<str_vec>( desc_and_help_of_option_name | map_keys );
}

/// \brief Construct a detail_help_options_block from a map from option name to a pair of description and help message
detail_help_options_block::detail_help_options_block(str_str_str_pair_map arg_desc_and_help_of_option_name ///< A map of option name (string) to a pair of description (string) and help message (string)
                                                     ) : desc_and_help_of_option_name { std::move( arg_desc_and_help_of_option_name ) } {
}

/// \brief Whether the options specified in this detail_help_options_block are requesting any help strings
bool detail_help_options_block::has_help_string() const {
	for (const str_bool_pair &value : values) {
		if (value.second) {
			return true;
		}
	}
	return false;
}

/// \brief Return the help string that the options specified in this detail_help_options_block are requesting
///
/// \pre has_help_string() should be true else an invalid_argument_exception will be thrown
///
/// \returns The help string that has been requested
string detail_help_options_block::help_string() const {
	for (const str_bool_pair &value : values) {
		if (value.second) {
			return desc_and_help_of_option_name.find(value.first)->second.second;;
		}
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to find any detail help options that were requested"));

	return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
}
