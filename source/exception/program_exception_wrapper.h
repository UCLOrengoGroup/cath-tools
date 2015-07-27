/// \file
/// \brief The program_exception_wrapper class header

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

#ifndef PROGRAM_EXCEPTION_WRAPPER_H_INCLUDED
#define PROGRAM_EXCEPTION_WRAPPER_H_INCLUDED

#include <boost/log/utility/setup/console.hpp>

#include "common/logger.h"

#include <iostream>

using sink_t    = boost::log::sinks::synchronous_sink<boost::log::sinks::basic_text_ostream_backend<char> > ;
using sink_sptr = boost::shared_ptr<sink_t>;

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
			/// \brief A shared_ptr to a Boost Log sink that is added to the core for the duration of the program run
			sink_sptr boost_log_sink_sptr;

			void output_catch_context(std::ostream &, const char * const) const;

			virtual std::string do_get_program_name() const = 0;
			virtual void do_run_program(int, char * []) = 0;

		public:
			program_exception_wrapper() = default;
			virtual ~program_exception_wrapper() noexcept;

			/// \brief Specify that the copy-ctor shouldn't be used
			program_exception_wrapper(const program_exception_wrapper &) = delete;
			/// \brief Specify that the move-ctor shouldn't be used
			program_exception_wrapper(program_exception_wrapper &&) = delete;
			/// \brief Specify that the copy-assign shouldn't be used
			program_exception_wrapper & operator=(const program_exception_wrapper &) = delete;
			/// \brief Specify that the move-assign shouldn't be used
			program_exception_wrapper & operator=(program_exception_wrapper &&) = delete;

			int run_program(int,
			                char * [],
			                std::ostream &arg_os = std::cerr);
			sink_sptr get_sink_ptr();
		};

	}
}

#endif
