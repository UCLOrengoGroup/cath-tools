/// \file
/// \brief The log_to_ostream_guard class definitions

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

#include <iostream>

#include <spdlog/sinks/ostream_sink.h>

#include "cath/common/boost_addenda/log/log_to_ostream_guard.hpp"

using namespace ::cath;

using ::std::make_shared;
using ::std::ostream;
using ::std::shared_ptr;

namespace sinks = ::spdlog::sinks;

/// \brief Ctor for log_to_ostream_guard
///
/// \param prm_ostream The ostream to which logging should temporarily be redirected for the log_to_ostream_guard's lifetime
log_to_ostream_guard::log_to_ostream_guard( ostream &prm_ostream ) : logger_shptr{ ::spdlog::default_logger() } {
	// The previous logger is being stored in logger_shptr
	// So now set a new logger to log to prm_ostream
	::spdlog::set_default_logger( make_shared<::spdlog::logger>(
	  "", shared_ptr<::sinks::sink>( make_shared<::sinks::ostream_sink_mt>( prm_ostream ) ) ) );

	// Change the format to just contain the message
	::spdlog::default_logger()->set_pattern( "%v" );
}

/// \brief Virtual empty dtor for log_to_ostream_guard
log_to_ostream_guard::~log_to_ostream_guard() noexcept {
	// If this log_to_ostream_guard has changed the default logger, reset it
	try {
		reset_default_logger();
	}
	// ...but don't let any exceptions escape the destructor
	catch ( ... ) {
	}
}

/// \brief Reset the default logger
void log_to_ostream_guard::reset_default_logger() {
	if ( logger_shptr ) {
		::spdlog::set_default_logger( logger_shptr );
		logger_shptr.reset();
	}
}
