/// \file
/// \brief The misc_help_version_options_block class definitions

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

#include "misc_help_version_options_block.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/optional.hpp>

#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/clone/make_uptr_clone.h"
#include "exception/invalid_argument_exception.h"

using namespace boost::algorithm;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::algorithm::is_any_of;
using boost::algorithm::join;
using boost::none;

/// The program's current version
const string misc_help_version_options_block::CATH_BINARIES_VERSION      ( "v0.12.4"            );

/// The program's most recent update date
const string misc_help_version_options_block::CATH_BINARIES_VERSION_DATE ( "19th January 2015"    );

/// \brief The option name for the help option
const string misc_help_version_options_block::PO_HELP   ( "help"    );

/// \brief The option name for the version option
const string misc_help_version_options_block::PO_VERSION( "version" );

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
void misc_help_version_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                            ) {
	arg_desc.add_options()
		((PO_HELP    + ",h").c_str(), bool_switch( &help    )->default_value( false ), "Output help message"        )
		((PO_VERSION + ",v").c_str(), bool_switch( &version )->default_value( false ), "Output version information" );
}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// This is a concrete definition of a virtual method that's pure in options_block
///
/// At present, this always accepts
opt_str misc_help_version_options_block::do_invalid_string() const {
	return none;
}

/// \brief Get the help flag
bool misc_help_version_options_block::get_help() const {
	return help;
}

/// \brief Get the version flag
bool misc_help_version_options_block::get_version() const {
	return version;
}

/// \brief Generate the help string
string misc_help_version_options_block::get_help_string(const options_description &arg_visible_program_options, ///< The full options_description of visible options
                                                        const string              &arg_help_message_prefix,     ///< The prefix to prepend to the output of the options_description
                                                        const string              &arg_help_message_suffix      ///< The suffix to append to the output of the options_description
                                                        ) const {
	ostringstream help_ss;
	help_ss << arg_help_message_prefix << endl;
	help_ss << arg_visible_program_options << endl;
	help_ss << arg_help_message_suffix;
	return help_ss.str();
}

/// \brief Generate the version string
string misc_help_version_options_block::get_version_string(const string &arg_program_name,       ///< The name of the program
                                                           const string &arg_program_description ///< A description of the program
                                                           ) const {
	// Construct a copy of the string that indents all of the lines
	const string  indent_string("   ");
	const str_vec desc_lines = split_build<str_vec>( arg_program_description, is_any_of( "\n" ) );
	const string  indented_description = indent_string + join( desc_lines, "\n" + indent_string );

	ostringstream version_ss;
	version_ss << "Overview"                                                                                 << endl;
	version_ss << "--------"                                                                                 << endl;
	version_ss << "   Name                  : " << arg_program_name                                          << endl;
	version_ss << "   Version               : " << CATH_BINARIES_VERSION                                     << endl;
	version_ss << "   Version date          : " << CATH_BINARIES_VERSION_DATE                                << endl;
	version_ss << endl;
	version_ss << "Build details"                                                                            << endl;
	version_ss << "-------------"                                                                            << endl;
	version_ss << "   Platform              : " << BOOST_PLATFORM                                            << endl;
	version_ss << "   Compiler              : " << BOOST_COMPILER                                            << endl;
	version_ss << "   Library               : " << BOOST_STDLIB                                              << endl;
	version_ss << "   Boost version         : " << BOOST_LIB_VERSION                                         << endl;
//	version_ss << "   Boost Program Options : " << BOOST_PROGRAM_OPTIONS_VERSION                             << endl;
#ifdef BUILD_BRANCH_NAME
	version_ss << "   Build branch          : " << BUILD_BRANCH_NAME                                         << endl;
#endif
#ifdef BUILD_REVISION_NUMBER
	version_ss << "   Build revision number : " << BUILD_REVISION_NUMBER                                     << endl;
#endif
	version_ss << "   Compile time          : " << __DATE__ << " " << __TIME__ << " (of " << __FILE__ << ")" << endl;
	version_ss << endl;
	version_ss << "Description"                                                                              << endl;
	version_ss << "-----------"                                                                              << endl;
	version_ss << indented_description                                                                       << endl;
	return version_ss.str();
}
