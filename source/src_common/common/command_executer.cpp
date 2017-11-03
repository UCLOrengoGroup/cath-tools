/// \file
/// \brief The command_executer class definitions

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#include "command_executer.hpp"

#include <boost/range/join.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/argc_argv_faker.hpp"
#include "common/exception/runtime_error_exception.hpp"

#include <sys/wait.h>

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::filesystem::path;
using boost::range::join;

/// \brief Execute a system command, allowing searching within PATH but not using the shell to process the command
///
/// This uses execvp, rather than system() because system()'s use of shell processing is very dangerous,
/// especially as the user has been allowed to specify the executable here. See:
///  https://www.securecoding.cert.org/confluence/display/cplusplus/ENV04-CPP.+Do+not+call+system%28%29+if+you+do+not+need+a+command+processor
///
/// \retval Success whether the call executed successfully
bool command_executer::execute(const path    &arg_command,  ///< TODOCUMENT
                               const str_vec &arg_arguments ///< TODOCUMENT
                               ) {
	// Prepare a vector of the command and the arguments
	const auto command_range = { arg_command.string() };
	const auto arguments     = copy_build<str_vec>( join( command_range, arg_arguments ) );
	argc_argv_faker arguments_faker( arguments );

	// Fork so that on half can become the command and the other half can wait for it and then continue
	const pid_t pid = fork();
	if (pid == -1) {
		BOOST_THROW_EXCEPTION(runtime_error_exception( "Unable to fork in command_executer::execute() to execute " + arg_command.string() ));
	}
	else if (pid == 0) {
		if ( execvp( arg_command.string().c_str(), arguments_faker.get_argv() ) == -1) {
			perror( ("Error executing " + arg_command.string()).c_str() );
			_exit(127);
		}
	}
	else {
		int status;
		pid_t ret;
		while ((ret = waitpid(pid, &status, 0)) == -1) {
			if (errno != EINTR) {
				BOOST_THROW_EXCEPTION(runtime_error_exception( "Error executing " + arg_command.string()));
				return false;
			}
		}
//		cerr << "Ret was    : " << ret    << endl;
//		cerr << "status was : " << status << endl;
		if (status != 0) {
			return false;
		}
	}
	return true;
}
