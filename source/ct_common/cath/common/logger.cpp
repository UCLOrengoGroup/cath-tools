/// \file
/// \brief The logger class definitions

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

#include "logger.hpp"

#include <boost/log/trivial.hpp>

#include <cstdlib>
#include <iostream>

using namespace ::cath;
using namespace ::std;

/// \brief TODOCUMENT
void logger::log_and_exit(const return_code     &prm_return_code,    ///< TODOCUMENT
                          const string          &exit_message,       ///< TODOCUMENT
                          const ostream_ref_opt &prm_message_ostream ///< TODOCUMENT
                          ) {
	if ( prm_message_ostream ) {
		prm_message_ostream->get() << exit_message;
	}
	else {
		BOOST_LOG_TRIVIAL( error ) << exit_message;
	}
	exit( static_cast<int>( prm_return_code ) );
}
