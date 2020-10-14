/// \file
/// \brief The program_exception_wrapper class definitions

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

#include "program_exception_wrapper.hpp"

#include <boost/core/demangle.hpp>
#include <boost/exception/all.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/numeric/conversion/cast.hpp>

using namespace ::cath::common;
using namespace ::std;

using ::boost::core::demangle;
using ::boost::log::expressions::format_date_time;
using ::boost::log::expressions::message;
using ::boost::log::expressions::smessage;
using ::boost::log::expressions::stream;
using ::boost::log::keywords::filter;
using ::boost::log::keywords::format;
using ::boost::log::trivial::severity;
using ::boost::log::trivial::warning;

/// \brief A simple private function to provide a standard way of outputting the context of a catch to a stream.
void program_exception_wrapper::output_catch_context(ostream            &prm_os,
                                                     const char * const  prm_program_name
                                                     ) const {
	prm_os << "Whilst running program " << prm_program_name;
	prm_os << " (via a program_exception_wrapper with typeid: \"" << demangle( typeid( *this ).name() ) << "\")";
}

/// \brief Virtual destructor for program_exception_wrapper.
program_exception_wrapper::~program_exception_wrapper() noexcept {
	// If this program_exception_wrapper has added a sink to the Boost Log core
	// but wasn't able to remove it (presumably because an exception was thrown),
	// then try to remove it now
	try {
		if ( boost_log_sink_sptr ) {
			boost::log::core::get()->remove_sink( boost_log_sink_sptr );
			boost_log_sink_sptr.reset();
		}
	}
	// ...but don't let any exceptions escape the destructor
	catch (...) {
	}
}

/// \brief A NVI wrapper that calls the virtual do_run_program() and catches and describes any exceptions
///
/// To create a new program, write a concrete class that derives from program_exception_wrapper and that
/// contains the program's implementation in the do_run_program() method.
///
/// Then just write a main that calls run_program on an instance of that class.

int program_exception_wrapper::run_program(int      prm_c,   ///< The main()-style argc parameter
                                           char *   prm_v[], ///< The main()-style argv parameter
                                           ostream &prm_os   ///< The ostream to which problems should be reported [Default: cerr (ie stderr)]
                                           ) {
	// Check that some arguments have been provided
	if (prm_c <= 0) {
		prm_os << "Cannot run_program() without a strictly positive number of arguments (the program name should be the first)" << endl;
		return static_cast<int>( logger::return_code::GENERIC_FAILURE_RETURN_CODE );
	}

	// Try calling do_run_program with the arguments.
	// Catch any errors and describe them to prm_os
	try {
		// If using Boost Log, then add a sink that writes to stderr, rather than using the default stdout sink
		boost_log_sink_sptr = boost::log::add_console_log(
			cerr,
			format = (
				stream          << format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S.%f")
//	                   << "] [" << attr<boost::log::attributes::current_thread_id::value_type >("ThreadID")
					   << " ["  << do_get_program_name()
					   << "|\033[1m"  << left << setw( 7 ) << setfill(' ') << severity
					   << "\033[0m] " << smessage
			),
			filter = ( severity >= warning )
		);
		boost::log::add_common_attributes();

		do_run_program(prm_c, prm_v);

		// If this has added a sink to the Boost Log core then try to remove it
		if ( boost_log_sink_sptr ) {
			boost::log::core::get()->remove_sink( boost_log_sink_sptr );
			boost_log_sink_sptr.reset();
		}
	}
	catch (const boost::exception &e) {
		output_catch_context(prm_os, prm_v[ 0 ] );
		prm_os << ", caught a boost::exception:\n";
		prm_os << diagnostic_information(e);
		prm_os << "\n" << endl;
		return static_cast<int>( logger::return_code::GENERIC_FAILURE_RETURN_CODE );
	}
	catch (const std::exception &e) {
		output_catch_context(prm_os, prm_v[ 0 ] );
		prm_os << ", caught a std::exception:\n";
		prm_os << e.what();
		prm_os << "\n" << endl;
		return static_cast<int>( logger::return_code::GENERIC_FAILURE_RETURN_CODE );
	}
	catch (...) {
		output_catch_context(prm_os, prm_v[ 0 ] );
		prm_os << ", caught an exception of unrecognised type";
		prm_os << "\n" << endl;
		return static_cast<int>( logger::return_code::GENERIC_FAILURE_RETURN_CODE );
	}

	return static_cast<int>( logger::return_code::SUCCESS );
}

sink_sptr program_exception_wrapper::get_sink_ptr() {
	return boost_log_sink_sptr;
}
