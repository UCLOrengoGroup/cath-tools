/// \file
/// \brief The program_exception_wrapper class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROGRAM_EXCEPTION_WRAPPER_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROGRAM_EXCEPTION_WRAPPER_HPP

#include <spdlog/spdlog.h>

#include <iostream>
#include <string>

namespace cath {
	namespace common {

		/// \brief This is an ABC from which programs can derive to get standard last-chance exception handling
		///
		/// If you're writing a new program, create a class that derives from this one
		/// and put the actions of the program inside the do_run_program() method.
		/// Then just make your main() function call run_program() on an instance of your class
		/// to get standard last-chance exception handling for free
		///
		/// Note that this design won't catch any exceptions that are thrown before do_run_program() is called.
		class program_exception_wrapper {
		  private:
			/// A holder for the default_logger that gets switched out during the lifetime of the run
			::std::shared_ptr<::spdlog::logger> logger_shptr;

			void output_catch_context( ::std::ostream &, const char * ) const;

			[[nodiscard]] virtual ::std::string do_get_program_name() const = 0;

			virtual void do_run_program( int, char *[] ) = 0;

			void reset_default_logger() noexcept;

		  public:
			program_exception_wrapper() = default;
			virtual ~program_exception_wrapper() noexcept;

			/// \brief Block the copy-ctor
			program_exception_wrapper( const program_exception_wrapper & ) = delete;
			/// \brief Block the move-ctor
			program_exception_wrapper( program_exception_wrapper && ) = delete;
			/// \brief Block the copy-assign
			program_exception_wrapper &operator=( const program_exception_wrapper & ) = delete;
			/// \brief Block the move-assign
			program_exception_wrapper &operator=( program_exception_wrapper && ) = delete;

			int run_program( int, char *[], ::std::ostream & = ::std::cerr );
		};

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROGRAM_EXCEPTION_WRAPPER_HPP
