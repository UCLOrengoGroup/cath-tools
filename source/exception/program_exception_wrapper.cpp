/// \file
/// \brief The program_exception_wrapper class definitions

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

#include "program_exception_wrapper.h"

#include <boost/exception/all.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/numeric/conversion/cast.hpp>

using namespace boost::log;
using namespace boost::log::expressions;
using namespace boost::log::trivial;
//using namespace cath;
using namespace cath::common;
using namespace std;

using boost::numeric_cast;

/// \brief A simple private function to provide a standard way of outputting the context of a catch to a stream.
void program_exception_wrapper::output_catch_context(ostream            &arg_os,
                                                     const char * const  arg_program_name
                                                     ) const {
	arg_os << "Whilst running program " << arg_program_name;
	arg_os << " (via a program_exception_wrapper with typeid: \"" << typeid( *this ).name() << "\")";
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

int program_exception_wrapper::run_program(int      arg_c,   ///< The main()-style argc parameter
                                           char *   arg_v[], ///< The main()-style argv parameter
                                           ostream &arg_os   ///< The ostream to which problems should be reported [Default: cerr (ie stderr)]
                                           ) {
	// Check that some arguments have been provided
	if (arg_c <= 0) {
		arg_os << "Cannot run_program() without a strictly positive number of arguments (the program name should be the first)" << endl;
		return static_cast<int>( logger::return_code::EXCEPTION_WITHOUT_SPECIFIC_RETCODE );
	}

	// Try calling do_run_program with the arguments.
	// Catch any errors and describe them to arg_os
	try {
		// If using Boost Log, then add a sink that writes to stderr, rather than using the default stdout sink
		boost_log_sink_sptr = boost::log::add_console_log(
			cerr,
			keywords::format = (
				stream          << format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S.%f")
//	                   << "] [" << attr<boost::log::attributes::current_thread_id::value_type >("ThreadID")
					   << " ["  << do_get_program_name()
					   << "|"   << left << setw( 7 ) << setfill(' ') << trivial::severity
					   << "] "  << smessage
			),
			keywords::filter = ( severity >= warning )
		);
		boost::log::add_common_attributes();

		do_run_program(arg_c, arg_v);

		// If this has added a sink to the Boost Log core then try to remove it
		if ( boost_log_sink_sptr ) {
			boost::log::core::get()->remove_sink( boost_log_sink_sptr );
			boost_log_sink_sptr.reset();
		}
	}
	catch (const boost::exception &e) {
		output_catch_context(arg_os, arg_v[ 0 ] );
		arg_os << ", caught a boost::exception:\n";
		arg_os << diagnostic_information(e);
		arg_os << "\n" << endl;
		return static_cast<int>( logger::return_code::EXCEPTION_WITHOUT_SPECIFIC_RETCODE );
	}
	catch (const std::exception &e) {
		output_catch_context(arg_os, arg_v[ 0 ] );
		arg_os << ", caught a std::exception:\n";
		arg_os << e.what();
		arg_os << "\n" << endl;
		return static_cast<int>( logger::return_code::EXCEPTION_WITHOUT_SPECIFIC_RETCODE );
	}
	catch (...) {
		output_catch_context(arg_os, arg_v[ 0 ] );
		arg_os << ", caught an exception of unrecognised type";
		arg_os << "\n" << endl;
		return static_cast<int>( logger::return_code::EXCEPTION_WITHOUT_SPECIFIC_RETCODE );
	}

	return static_cast<int>( logger::return_code::SUCCESS );
}

sink_sptr program_exception_wrapper::get_sink_ptr() {
	return boost_log_sink_sptr;
}
