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

#include "log_to_ostream_guard.hpp"

#include <boost/core/null_deleter.hpp>
#include <boost/shared_ptr.hpp>

using namespace cath;
using namespace std;

using boost::null_deleter;

/// \brief Ctor for log_to_ostream_guard
log_to_ostream_guard::log_to_ostream_guard(ostream &arg_ostream ///< TODOCUMENT
                                           ) {
	// Construct a sink
    boost_log_sink_bsptr = boost::make_shared<sink_t>();

    // Get a (non-deleting) shared_ptr to the ostream and add it to the sink
    boost::shared_ptr<ostream> stream_bsptr( &arg_ostream, null_deleter() );
    boost_log_sink_bsptr->locked_backend()->add_stream( stream_bsptr );

    // Register the sink in the logging core
    boost::log::core::get()->add_sink( boost_log_sink_bsptr );
}

/// \brief Virtual empty dtor for log_to_ostream_guard
log_to_ostream_guard::~log_to_ostream_guard() noexcept {
	// If this log_to_ostream_guard has added a sink to the Boost Log core
	// then try to remove it now
	try {
		remove_log_sink();
	}
	// ...but don't let any exceptions escape the destructor
	catch (...) {
	}
}

/// \brief TODOCUMENT
void log_to_ostream_guard::remove_log_sink() {
	if ( boost_log_sink_bsptr ) {
		boost::log::core::get()->remove_sink( boost_log_sink_bsptr );
		boost_log_sink_bsptr.reset();
	}
}
