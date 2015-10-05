/// \file
/// \brief The executable_options class definitions

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

#include "executable_options.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/optional.hpp>

#include "common/argc_argv_faker.h"
#include "common/file/open_fstream.h"
#include "exception/invalid_argument_exception.h"
#include "options/executable/env_var_option_name_handler.h"
#include "options/options_block/data_dirs_options_block.h"
#include "options/options_block/misc_help_version_options_block.h"

#include <fstream>
#include <iostream>
#include <iosfwd>

using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

constexpr size_t executable_options::DEFAULT_PROG_OPS_LINE_LENGTH;

/// \brief The prefix for "global" environment variables to be respected by all executables using this code.
///
/// These "global" environment variables are overridden by command line options but override the configuration file.
const string       executable_options::CATH_TOOLS_ENVIRONMENT_VARIABLE_PREFIX("CATH_TOOLS_");

/// \brief The name of the "global" configuration file to be respected by all executables using this code.
///
/// This "global" file is overridden by environment variables or command line options
const path         executable_options::CATH_TOOLS_CONF_FILE                  ("cath-tools.conf");

/// \brief The path through which to search for the "global" configuration file.
const path_vec executable_options::CATH_TOOLS_CONF_FILE_SEARCH_PATH = { ".", "~/.cath" };

/// \brief Default implementation of get_positional_options() that doesn't set any positional options
///
/// This is a virtual function that may be overridden by any concrete, derived classes
positional_options_description executable_options::get_positional_options() {
	return positional_options_description();
}

/// \brief Return a standard string that can be appended to any options error message
///
/// The current string is of the form:
///
/// > Try 'cath-ssap --help' for usage information.
string executable_options::get_standard_usage_error_string() const {
	return "Try '"
	       + get_program_name()
	       + " --"
	       + misc_help_version_options_block::PO_HELP
	       + "' for usage information.";
}

/// \brief An NVI pass-though method to get the name of the program from do_get_program_name()
string executable_options::get_program_name() const {
	return do_get_program_name();
}

/// \brief Get the Boost variables map that has been populated by parsing options
///
/// \pre The executable_options should have already parsed options,
///      else invalid_argument_exception will be thrown
///
/// Many executables won't need to use this method but it is needed to query
/// things like whether options had values specified or were defaulted() etc
const variables_map & executable_options::get_variables_map() const {
	// Check the options have been processed
	if (!processed_options) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get variables_map because the options haven't yet been processed"));
	}
	return vm;
}

/// \brief Add an options block to the list to be handled by this executable_options
///
/// \pre The executable_options should not yet have parsed options,
///      else invalid_argument_exception will be thrown
///
/// Note that this registers the options_block but doesn't take ownership of it
/// The client should
void executable_options::add_options_block(options_block &arg_options_block ///< The options block to register
                                           ) {
	if (processed_options) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot add an options_block once the options have been processed"));
	}
	all_options_blocks.push_back(&arg_options_block);
}

/// \brief Attempt to parse the specified options
///
/// \pre All previously added options_blocks must still exist
///
/// \pre The two parameters must be in the standard argc/argv configuration
///
/// Note that the two input parameters are taken const
/// (including the argv being (a pointer to) an array of const pointers *to const*)
/// so they will remain untouched.
///
/// \post The options will be parsed and ready for querying
void executable_options::parse_options(const int          &argc,  ///< The argc from command line parameters
                                       const char * const  argv[] ///< The argv from command line parameters
                                       ) {
	// Check the options haven't already been processed
	if (processed_options) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot process options once they have already been processed"));
	}

	// Create two options_description, one complete and another containing all visible options
	options_description full_po_desc   ( DEFAULT_PROG_OPS_LINE_LENGTH );
	options_description visible_po_desc( DEFAULT_PROG_OPS_LINE_LENGTH );

	// Frustratingly, Boost 1.41 (as used by orengobuild64) won't accept a const argv, so
	// it's necessary to construct a fake argc/argv here
	argc_argv_faker fake_argc_argv(argc, argv);
	int new_argc      = fake_argc_argv.get_argc();
	char * * new_argv = fake_argc_argv.get_argv();

	// Add each of the options_blocks to po_desc
	for (options_block * const options_block_ptr : all_options_blocks) {
		const options_description hidden_opts  = options_block_ptr->get_hidden_options_description(  DEFAULT_PROG_OPS_LINE_LENGTH );
		const options_description visible_opts = options_block_ptr->get_visible_options_description( DEFAULT_PROG_OPS_LINE_LENGTH );
		full_po_desc.add(    hidden_opts  );
		full_po_desc.add(    visible_opts );
		visible_po_desc.add( visible_opts );
	}

	// Attempt the parses and catch any exceptions
	//
	// The parses are performed in decreasing order of precedence
	// (ie options specified via the command line should take precedence over those
	//  specified via environment variables so it comes first)
	//
	// \todo If this gets any more complicated then consider putting each of these
	//       different parses into different objects (presumably of classes deriving from
	//       a single ABC).
	string parsing_approach = "";
	try {
		// Parse options from the command line
		parsing_approach = "from the command line";
		const positional_options_description positionals = get_positional_options();
		processed_options = true;
		store(
			command_line_parser(
				new_argc,
				new_argv
			).options(
				full_po_desc
			).positional(
				positionals
			).run(),
			vm
		);

		// Parse any environment variables prefixed with "CATH_TOOLS_"
		// and just silently ignore any unknown options
		parsing_approach = "from global environment variables with prefix " + CATH_TOOLS_ENVIRONMENT_VARIABLE_PREFIX;
		//! [Using env_var_option_name_handler]
		store(
			parse_environment(
				full_po_desc,
				env_var_option_name_handler(
					CATH_TOOLS_ENVIRONMENT_VARIABLE_PREFIX,
					true,
					full_po_desc
				)
			),
			vm
		);
		//! [Using env_var_option_name_handler]

		// Parse any configuration file called cath_tools.conf
		const path located_cath_tools_conf_file = find_file(CATH_TOOLS_CONF_FILE_SEARCH_PATH, CATH_TOOLS_CONF_FILE.string());
		if (!located_cath_tools_conf_file.empty()) {
//			cerr << "Parsing configuration from file " << CATH_TOOLS_CONF_FILE << endl;
			parsing_approach = "from the global configuration file " + located_cath_tools_conf_file.string();
			ifstream config_file_stream;
			open_ifstream(config_file_stream, CATH_TOOLS_CONF_FILE);
			store(parse_config_file(config_file_stream, full_po_desc, true), vm);
			config_file_stream.close();
		}

		// All parsing is complete so call notify, which will trigger any
		// post-parsing hooks to get called
		notify(vm);
	}
	catch (std::exception &e) {
		error_or_help_string = get_program_name() + ": Error in parsing program options (" + parsing_approach + "): " + e.what();
	}
	catch (...) {
		error_or_help_string = get_program_name() + ": Caught an unrecognised exception whilst parsing program options (" + parsing_approach + ").\n";
	}

	// If no error has yet been encountered
	// and if any of the blocks return non-empty invalid_string() results,
	// then set error_or_help_string to the first
	if (error_or_help_string.empty()) {
		for (options_block * const options_block_ptr : all_options_blocks) {
			const opt_str block_invalid_string = options_block_ptr->invalid_string();
			if ( block_invalid_string ) {
				error_or_help_string = *block_invalid_string;
				break;
			}
		}
	}

	// If an error string has already arisen from the parsing or from an options block being invalid...
	if (!error_or_help_string.empty()) {
		// ...and if doesn't already end with the standard usage error string
		// (eg "Try 'cath-ssap --help' for usage information.")
		// then append it now.
		if (!ends_with(error_or_help_string, get_standard_usage_error_string())) {
			error_or_help_string += ("\n" + get_standard_usage_error_string());
		}
	}
	// Else there is no error string yet, so now let the derived class consider then
	// set it to whatever the derived class decides
	else {
		error_or_help_string = do_update_error_or_help_string( visible_po_desc );
	}
}

/// \brief Get any error or help string that was generated during parse
///
/// \pre parse_options() must have been successfully called on this object,
///      else this will throw an invalid_argument_exception
///
/// \returns Any error or help string generated during parsing,
///          or otherwise an empty string
string executable_options::get_error_or_help_string() const {
	if (!processed_options) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get error/help string because the options haven't yet been processed"));
	}
	return error_or_help_string;
}


