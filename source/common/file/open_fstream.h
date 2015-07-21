/// \file
/// \brief The open_ifstream / open_ofstream header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Binaries project and then tweaked, eg namespaced in cath)
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

#ifndef OPEN_FSTREAM_H_INCLUDED
#define OPEN_FSTREAM_H_INCLUDED

#include <boost/filesystem.hpp>

#include "exception/runtime_error_exception.h"

#include <iosfwd>

namespace cath {
	namespace common {
		namespace detail {

			enum class fstream_type {
				READING,
				WRITING
			};

			/// \brief Function used to implement open_ifstream and open_ofstream
			///
			/// This:
			/// - sets the exceptions to throw on failbit or badbit
			/// - tries to open the file
			/// - catches any failure and rethrows as a runtime_error_exception
			template <typename fstream_t>
			void open_fstream_impl(fstream_t                     &arg_fstream,     ///< TODOCUMENT
			                       const boost::filesystem::path &arg_filename,    ///< TODOCUMENT
			                       const std::ios_base::openmode &arg_mode,        ///< TODOCUMENT
			                       const fstream_type            &arg_fstream_type ///< TODOCUMENT
			                       ) {
				const std::string reading_or_writing_str = (arg_fstream_type == fstream_type::READING) ? "reading" : "writing";
				arg_fstream.exceptions(std::ios::failbit | std::ios::badbit);
				try {
					arg_fstream.open(arg_filename.string().c_str(), arg_mode);
				}
				// Catch any I/O exceptions
				catch (const std::exception &ex) {
			//		const std::string reading_or_writing_str = (arg_fstream_type == FILE_READING) ? "reading" : "writing";
					const std::string error_message(
						"Cannot open file \""
						+ arg_filename.string()
						+ "\" for "
						+ reading_or_writing_str
						+ " ["
						+ ex.what()
						+ "] "
					);
					perror(error_message.c_str());
					BOOST_THROW_EXCEPTION(cath::common::runtime_error_exception(error_message));
				};

				assert(arg_fstream.is_open());
				assert(arg_fstream.good());
			}
		}

		void open_ifstream(std::ifstream &,
		                   const boost::filesystem::path &,
		                   const std::ios_base::openmode &arg_mode = std::ios_base::in);

		void open_ofstream(std::ofstream &,
		                   const boost::filesystem::path &,
		                   const std::ios_base::openmode &arg_mode = std::ios_base::out);

	}
}

#endif
