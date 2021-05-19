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

#include <boost/core/demangle.hpp>
#include <boost/exception/all.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "cath/common/logger.hpp"
#include "cath/common/program_exception_wrapper.hpp"

using namespace ::cath::common;
using namespace ::std;

using ::boost::core::demangle;

/// \brief A simple private function to provide a standard way of outputting the context of a catch to a stream.
void program_exception_wrapper::output_catch_context( ostream &prm_os, const char *const prm_program_name ) const {
	prm_os << "Whilst running program " << prm_program_name;
	prm_os << " (via a program_exception_wrapper with typeid: \"" << demangle( typeid( *this ).name() ) << "\")";
}

/// \brief TODOCUMENT
void program_exception_wrapper::reset_default_logger() noexcept {
	// If this program_exception_wrapper has stored a logger in logger_shptr,
	// return it as the default logger
	try {
		if ( logger_shptr ) {
			::spdlog::set_default_logger( logger_shptr );
			logger_shptr.reset();
		}
	}
	// ...but don't let any exceptions escape
	catch ( ... ) {
	}
}

/// \brief Virtual destructor for program_exception_wrapper.
program_exception_wrapper::~program_exception_wrapper() noexcept {
	reset_default_logger();
}

/// \brief A NVI wrapper that calls the virtual do_run_program() and catches and describes any exceptions
///
/// To create a new program, write a concrete class that derives from program_exception_wrapper and that
/// contains the program's implementation in the do_run_program() method.
///
/// Then just write a main that calls run_program on an instance of that class.
///
/// \param prm_c   The main()-style argc parameter
/// \param prm_v   The main()-style argv parameter
/// \param prm_os  The ostream to which problems should be reported [Default: cerr (ie stderr)]
int program_exception_wrapper::run_program( int prm_c, char *prm_v[], ostream &prm_os ) {
	// Check that some arguments have been provided
	if ( prm_c <= 0 ) {
		prm_os << "Cannot run_program() without a strictly positive number of arguments (the program name should be "
		          "the first)"
		       << endl;
		return static_cast<int>( logger::return_code::GENERIC_FAILURE_RETURN_CODE );
	}

	// Try calling do_run_program with the arguments.
	// Catch any errors and describe them to prm_os
	try {
		logger_shptr = ::spdlog::default_logger();
		// Add a sink that writes to stderr, rather than using the default stdout sink
		::spdlog::set_default_logger( ::spdlog::stderr_color_st( "program_exception_wrapper_logger" ) );
		::spdlog::default_logger()->set_level( spdlog::level::warn );
		::spdlog::default_logger()->set_pattern( ::fmt::format( "[%Y-%m-%d %H:%M:%S.%f] [%P:%t] [{}|%^%-8l%$] %v", do_get_program_name() ) );
		do_run_program( prm_c, prm_v );
	} catch ( const boost::exception &e ) {
		output_catch_context( prm_os, prm_v[ 0 ] );
		prm_os << ", caught a boost::exception:\n" << diagnostic_information( e ) << "\n" << endl;
		reset_default_logger();
		return static_cast<int>( logger::return_code::GENERIC_FAILURE_RETURN_CODE );
	} catch ( const std::exception &e ) {
		output_catch_context( prm_os, prm_v[ 0 ] );
		prm_os << ", caught a std::exception:\n" << e.what() << "\n" << endl;
		reset_default_logger();
		return static_cast<int>( logger::return_code::GENERIC_FAILURE_RETURN_CODE );
	} catch ( ... ) {
		output_catch_context( prm_os, prm_v[ 0 ] );
		prm_os << ", caught an exception of unrecognised type\n" << endl;
		reset_default_logger();
		return static_cast<int>( logger::return_code::GENERIC_FAILURE_RETURN_CODE );
	}
	reset_default_logger();
	return static_cast<int>( logger::return_code::SUCCESS );
}
