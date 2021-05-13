/// \file
/// \brief The misc_help_version_options_block class definitions

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

#include "misc_help_version_options_block.hpp"

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/external_info/cath_tools_git_version.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

using ::boost::lexical_cast;
using ::boost::program_options::bool_switch;
using ::boost::program_options::options_description;
using ::boost::program_options::variables_map;
using ::std::nullopt;

/// \brief The option name for the hidden help option
const string misc_help_version_options_block::PO_HIDDEN_HELP          { "hidden-help"            };

/// \brief The option name for the help option
const string misc_help_version_options_block::PO_HELP                 { "help"                   };

/// \brief The option name for the version option
const string misc_help_version_options_block::PO_VERSION              { "version"                };

/// \brief The single-character for the help option
constexpr char misc_help_version_options_block::PO_CHAR_HELP;

/// \brief The single-character for the version option
constexpr char misc_help_version_options_block::PO_CHAR_VERSION;

/// \brief A standard do_clone method
///
/// This is a concrete definition of a virtual method that's pure in options_block
unique_ptr<options_block> misc_help_version_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
///
/// This is a concrete definition of a virtual method that's pure in options_block
string misc_help_version_options_block::do_get_block_name() const {
	return "Miscellaneous";
}

/// \brief Add this block's options to the provided options_description
///
/// This is a concrete definition of a virtual method that's pure in options_block
void misc_help_version_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                            const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                            ) {
	const string PO_HELP_W_CHAR    = PO_HELP    + ',' + PO_CHAR_HELP;
	const string PO_VERSION_W_CHAR = PO_VERSION + ',' + PO_CHAR_VERSION;
	prm_desc.add_options()
		( PO_HELP_W_CHAR    . c_str(), bool_switch( &help    )->default_value( false ), "Output help message"        )
		( PO_VERSION_W_CHAR . c_str(), bool_switch( &version )->default_value( false ), "Output version information" );
}

/// \brief Add a hidden option to the options_description for the hidden help option
void misc_help_version_options_block::do_add_hidden_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                           const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                           ) {
	prm_desc.add_options()
		( ( PO_HIDDEN_HELP ).c_str(), bool_switch( &hidden_help )->default_value( false ), "Output help message including all hidden options" );
}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// This is a concrete definition of a virtual method that's pure in options_block
///
/// At present, this always accepts
str_opt misc_help_version_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                           ) const {
	return nullopt;
}

/// \brief Return all options names for this block
str_vec misc_help_version_options_block::do_get_all_options_names() const {
	return {
		misc_help_version_options_block::PO_HELP,
		misc_help_version_options_block::PO_VERSION,
	};
}

/// \brief Get the hidden_help flag
const bool & misc_help_version_options_block::get_hidden_help() const {
	return hidden_help;
}

/// \brief Get the help flag
const bool & misc_help_version_options_block::get_help() const {
	return help;
}

/// \brief Get the version flag
const bool & misc_help_version_options_block::get_version() const {
	return version;
}

/// \brief Generate the help string
string misc_help_version_options_block::get_help_string(const options_description &prm_visible_program_options, ///< The full options_description of visible options
                                                        const string              &prm_help_message_prefix,     ///< The prefix to prepend to the output of the options_description
                                                        const string              &prm_help_message_suffix      ///< The suffix to append to the output of the options_description
                                                        ) {
	return prm_help_message_prefix + "\n"
		+ lexical_cast<string>( prm_visible_program_options )
		+ prm_help_message_suffix;
}

// Since Clang and GCC indicate the address sanitizer (ASAN) in different ways,
// set the GCC flag if ASAN's detected under Clang
#ifdef __clang__
#if __has_feature(address_sanitizer)
#define __SANITIZE_ADDRESS__
#endif
#endif

/// \brief Generate the version string
string misc_help_version_options_block::get_version_string(const string &prm_program_name,       ///< The name of the program
                                                           const string &prm_program_description ///< A description of the program
                                                           ) {
	return "============\n"
		+ prm_program_name + " " + cath_tools_git_version() + " [" + cath_tools_git_date() + "]\n"
		+ "============\n"
		+ "\n"
		+ prm_program_description + "\n"
		+ "\n"
		+ "Build\n"
		+ "-----\n"
		+ "   "       + __DATE__ + " " + __TIME__  + "\n"
		+ "   "       + BOOST_COMPILER
#if defined(__SANITIZE_ADDRESS__)
		                               + " [ASAN]"
#endif
		                                           + "\n"
		+ "   "       + BOOST_STDLIB               + "\n"
		+ "   Boost " + BOOST_LIB_VERSION          + "\n";
}
