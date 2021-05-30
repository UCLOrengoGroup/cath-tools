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

#include "executable_options.hpp"

#include <filesystem>
#include <fstream>
#include <iosfwd>
#include <iostream>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <fmt/core.h>

#include "cath/common/argc_argv_faker.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/file/find_file.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/optional/make_optional_if.hpp"
#include "cath/common/test_or_exe_run_mode.hpp"
#include "cath/options/executable/env_var_option_name_handler.hpp"
#include "cath/options/options_block/misc_help_version_options_block.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;

using ::boost::program_options::command_line_parser;
using ::boost::program_options::options_description;
using ::boost::program_options::positional_options_description;
using ::boost::program_options::variables_map;
using ::boost::trim_right_copy;
using ::std::filesystem::path;
using ::std::ifstream;
using ::std::literals::string_literals::operator""s;
using ::std::make_optional;
using ::std::string;
using ::std::string_view;

/// \brief The name of the "global" configuration file to be respected by all executables using this code.
///
/// This "global" file is overridden by environment variables or command line options
static const path CATH_TOOLS_CONF_FILE( "cath-tools.conf" );

/// \brief The path through which to search for the "global" configuration file.
path_vec executable_options::CATH_TOOLS_CONF_FILE_SEARCH_PATH () {
	return path_vec{ ".", "~/.cath" };
}

/// \brief Default implementation of get_positional_options() that doesn't set any positional options
///
/// This is a virtual function that may be overridden by any concrete, derived classes
positional_options_description executable_options::get_positional_options() {
	return {};
}

/// \brief Return a standard string that can be appended to any options error message
///
/// The current string is of the form:
///
/// > Try 'cath-ssap --help' for usage information.
string executable_options::get_standard_usage_error_string() const {
	return ::fmt::format( "See '{} -{}' for usage.", get_program_name(), misc_help_version_options_block::PO_CHAR_HELP );
}

/// \brief An NVI pass-though method to get the name of the program from do_get_program_name()
string_view executable_options::get_program_name() const {
	return do_get_program_name();
}

/// \brief Get a string to prepend to the standard help
string executable_options::get_help_prefix_string() const {
	return do_get_help_prefix_string();
}

/// \brief Get a string to append to the standard help
string executable_options::get_help_suffix_string() const {
	return do_get_help_suffix_string()
		+ "\nPlease tell us your cath-tools bugs/suggestions : "
		+ "https://github.com/UCLOrengoGroup/cath-tools/issues/new";
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string executable_options::get_overview_string() const {
	return do_get_overview_string();
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
	if ( ! processed_options ) {
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
void executable_options::add_options_block(options_block &prm_options_block ///< The options block to register
                                           ) {
	if (processed_options) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot add an options_block once the options have been processed"));
	}
	all_options_blocks.emplace_back( prm_options_block );
}

/// \brief Add a string into the options usage
///
/// This is implemented by adding a no-options string_options_block.
/// The slight annoyance is that it will automatically append a ':'.
void executable_options::add_string(string prm_string ///< The string to add into the options usage
                                    ) {
	string_obj_blocks.emplace_back( ::std::move( prm_string ) );
	all_options_blocks.emplace_back( string_obj_blocks.back() );
}

/// \brief Add all the options of the specified options_block to the specified options_description
void executable_options::add_all_options_to_description(options_description &prm_options_description, ///< The options_description to which the options should be added
                                                        options_block       &prm_options_block,       ///< The options_block from which the options to be added should be drawn
                                                        const size_t        &prm_prog_ops_line_length ///< The length of the line to use
                                                        ) {
	prm_options_description.add( prm_options_block.get_visible_options_description( prm_prog_ops_line_length ) );
	const auto hidden_opts_block = prm_options_block.get_hidden_options_description ( prm_prog_ops_line_length );
	if ( ! hidden_opts_block.options().empty() ) {
		prm_options_description.add( hidden_opts_block );
	}
}

/// \brief Add the visible options of the specified options_block to the specified options_description
void executable_options::add_visble_options_to_description(options_description &prm_options_description, ///< The options_description to which the options should be added
                                                           options_block       &prm_options_block,       ///< The options_block from which the options to be added should be drawn
                                                           const size_t        &prm_prog_ops_line_length ///< The length of the line to use
                                                           ) {
	prm_options_description.add( prm_options_block.get_visible_options_description( prm_prog_ops_line_length ) );
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
void executable_options::parse_options(const int           &argc,             ///< The argc from command line parameters
                                       const char * const   argv[],           ///< The argv from command line parameters
                                       const parse_sources &prm_parse_sources ///< The sources from which options should be parsed
                                       ) {
	// Check the options haven't already been processed
	if ( processed_options ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot process options once they have already been processed"));
	}
	processed_options = true;

	// Create two options_description, one complete and another containing all visible options
	options_description full_po_desc   ( DEFAULT_PROG_OPS_LINE_LENGTH );
	options_description visible_po_desc( DEFAULT_PROG_OPS_LINE_LENGTH );

	// Frustratingly, Boost 1.41 (as used by orengobuild64) won't accept a const argv, so
	// it's necessary to construct a fake argc/argv here
	argc_argv_faker fake_argc_argv(argc, argv);
	int new_argc      = fake_argc_argv.get_argc();
	char * * new_argv = fake_argc_argv.get_argv();

	misc_help_version_options_block the_help_block;
	add_all_options_to_description   ( full_po_desc,    the_help_block, DEFAULT_PROG_OPS_LINE_LENGTH );
	add_visble_options_to_description( visible_po_desc, the_help_block, DEFAULT_PROG_OPS_LINE_LENGTH );

	// std::cerr << "full_po_desc    is " << full_po_desc    << "\n\n";

	// std::cerr << "visible_po_desc is " << visible_po_desc << "\n\n";

	// Attempt the standard parses and catch any exceptions
	//
	// \todo Consider putting each parse into different objects
	//       (presumably of classes deriving from a single ABC).

	// const auto parsed = try_parse(
	// 	command_line_parser( ac, av )
	// 		.options   ( action_od )
	// 		.positional( pod       )
	// 		.allow_unregistered(),
	// 	vm
	// );

	// return prog_opts_try(
	// 	[&] {
	// 		const auto parsed = prm_parser.run();
	// 		store( parsed, prm_vm );
	// 		notify( prm_vm );
	// 		return parsed;
	// 	}
	// );

	prog_opts_try(
		error_or_help_string,
		[&] {
			const auto parsed = command_line_parser( new_argc, new_argv )
				.options( full_po_desc )
				.allow_unregistered()
				.run();
			store( parsed, vm );
			notify( vm );
			// return parsed;
		},
		"[whilst parsing initial help/version options]"s
	);

	// std::cerr << "Done initial parse\n";

	// Add each of the options_blocks to po_desc
	for (const auto &the_options_block_ref : all_options_blocks) {
		add_all_options_to_description   ( full_po_desc,    the_options_block_ref.get(), DEFAULT_PROG_OPS_LINE_LENGTH );
		add_visble_options_to_description( visible_po_desc, the_options_block_ref.get(), DEFAULT_PROG_OPS_LINE_LENGTH );
	}

	// If (hidden) help was requested, then provide it
	if ( the_help_block.get_help() || the_help_block.get_hidden_help() ) {
		error_or_help_string = misc_help_version_options_block::get_help_string(
			( the_help_block.get_hidden_help() ? full_po_desc : visible_po_desc ),
			get_help_prefix_string(),
			get_help_suffix_string()
		);
		return;
	}

	// If version information was requested, then provide it
	if ( the_help_block.get_version() ) {
		// std::cerr << "Getting version\n";
		error_or_help_string =
		  misc_help_version_options_block::get_version_string( string( get_program_name() ), get_overview_string() );
		return;
	}

	// std::cerr << "On to main parses...\n";

	const positional_options_description positionals = get_positional_options();

	

	// Attempt the command line parse
	//
	// The remaining parses are performed in decreasing order of precedence
	// (ie options specified via the command line should take precedence over those
	//  specified via environment variables so it comes first)
	prog_opts_try(
		error_or_help_string,
		[&] {
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
		}
	);

	if ( prm_parse_sources == parse_sources::CMND_ENV_AND_FILE ) {

		// If running in test mode, complain that shouldn't be parsing env-vars or config files
		//
		// Note an alternative way to implement this is to get the global fixture to
		// detect and remove any CATH_TOOLS environment variables (though that doesn't
		// handle config files).
		if ( run_mode_flag::value == run_mode::TEST ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Requested to parse command line options from environment variables and configuration file but this is a bad idea because currently running in test mode"));
		}

		// Parse any environment variables prefixed with "CATH_TOOLS_"
		// and just silently ignore any unknown options
		//
		// The remaining parses are performed in decreasing order of precedence
		// (ie options specified via the command line should take precedence over those
		//  specified via environment variables so it comes first)
		prog_opts_try(
			error_or_help_string,
			[&] {
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
			},
			::fmt::format( "[whilst parsing from global environment variables with prefix {}]", CATH_TOOLS_ENVIRONMENT_VARIABLE_PREFIX )
		);

		// Parse any configuration file called cath_tools.conf
		const path located_cath_tools_conf_file = find_file( CATH_TOOLS_CONF_FILE_SEARCH_PATH(), CATH_TOOLS_CONF_FILE.string() );
		if ( ! located_cath_tools_conf_file.empty() ) {
//			cerr << "Parsing configuration from file " << CATH_TOOLS_CONF_FILE << endl;
			ifstream config_file_stream;
			prog_opts_try(
				error_or_help_string,
				[&] {
					open_ifstream( config_file_stream, CATH_TOOLS_CONF_FILE );
					store( parse_config_file( config_file_stream, full_po_desc, true ), vm );
					config_file_stream.close();
				},
				"[whilst parsing from the global configuration file " + located_cath_tools_conf_file.string() + "]"
			);
		}
	}

	// All parsing is complete so call notify, which will trigger any
	// post-parsing hooks to get called
	prog_opts_try(
		error_or_help_string,
		[&] { notify( vm ); }
	);

	if ( error_or_help_string ) {
		return;
	}

	// If no error has yet been encountered
	// and if any of the blocks return non-empty invalid_string() results,
	// then set error_or_help_string to the first
	for (const auto &options_block_ref : all_options_blocks) {
		const auto block_invalid_string = options_block_ref.get().invalid_string( vm );
		if ( block_invalid_string ) {
			prog_opts_try( error_or_help_string, [&] {
				BOOST_THROW_EXCEPTION(invalid_argument_exception(*block_invalid_string));
			} );
			return;
		}
	}

	// If there is still no error/help string yet, let the derived class have a chance to set it
	const auto concrete_error_or_help_string = do_get_error_or_help_string();
	if ( concrete_error_or_help_string ) {
		error_or_help_string = *concrete_error_or_help_string
			+ ( concrete_error_or_help_string->empty() ? ""s : "\n"s )
			+ get_standard_usage_error_string();
	}
}

/// \brief Get any error or help string that was generated during parse
///
/// \pre parse_options() must have been successfully called on this object,
///      else this will throw an invalid_argument_exception
///
/// \returns Any error or help string generated during parsing,
///          or otherwise an empty string
str_opt executable_options::get_error_or_help_string() const {
	if ( ! processed_options) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get error/help string because the options haven't yet been processed"));
	}
	return if_then_optional(
		static_cast<bool>( error_or_help_string ),
		trim_right_copy( *error_or_help_string ) + "\n"
	);
}


