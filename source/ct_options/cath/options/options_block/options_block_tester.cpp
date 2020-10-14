/// \file
/// \brief The options_block_tester class definitions

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

#include "options_block_tester.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>

#include "cath/common/argc_argv_faker.hpp"
#include "cath/options/options_block/options_block.hpp"

using namespace cath;
using namespace cath::opts;
using namespace std;

using boost::program_options::command_line_parser;
using boost::program_options::options_description;
using boost::program_options::variables_map;

const string options_block_tester::UNKNOWN_OPT  ( "--it-does-not-know-me"            );
const string options_block_tester::TEST_OPTION_1( "test_help_option_1"               );
const string options_block_tester::TEST_OPTION_2( "test_help_option_2"               );
const string options_block_tester::TEST_HELP_1  ( "This is the first piece of help"  );
const string options_block_tester::TEST_HELP_2  ( "This is the second piece of help" );

/// \brief Simple implementation function to return the specified list of parameters but with a dummy program name prepended
str_vec options_block_tester::prepend_dummy_program_name_copy(const str_vec &prm_params
                                                              ) {
	str_vec new_params( 1, "dummy_program_name" );
	new_params.reserve( prm_params.size() + 1 );
	for (const string &param : prm_params) {
		new_params.push_back( param );
	}
	return new_params;
}

/// \brief A simple method of parsing vector of option strings into an options_block for testing purposes
///
/// The options argument should just contain the parameters themselves;
/// this code will prepend a dummy program name to the beginning of them before use.
///
/// \todo Consider adding optional testing of positional parameters
void options_block_tester::parse_into_options_block(options_block &the_options_block, ///< The options_block to have options parsed into it
                                                    const str_vec &prm_options        ///< A vector of options strings to parse into the options_block (without the program name at the start - a dummy program name will be added)
	                                                ) {
	const str_vec opts_with_dummy_progname = prepend_dummy_program_name_copy( prm_options );
	options_description  po_desc           = the_options_block.get_all_options_description( 100 );
	argc_argv_faker      faked_arguments( opts_with_dummy_progname );

	// cerr << "Parsing from options : " << faked_arguments << endl;

	variables_map vm;
	store(
		command_line_parser(
			faked_arguments.get_argc(),
			faked_arguments.get_argv()
		).options(
			po_desc
		// ).positional(
		// 	positionals
		).run(),
		vm
	);

	// All parsing is complete so call notify, which will trigger any post-parsing hooks to get called
	notify( vm );
}

/// \brief TODOCUMENT
str_str_str_pair_map options_block_tester::TEST_DESC_AND_HELP_OF_OPTION_NAME() {
	return {
		{ TEST_OPTION_1, make_pair( "Description of " + TEST_OPTION_1, TEST_HELP_1 ), },
		{ TEST_OPTION_2, make_pair( "Description of " + TEST_OPTION_2, TEST_HELP_2 ), },
	};
}
