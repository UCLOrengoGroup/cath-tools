/// \file
/// \brief The open_ifstream / open_ofstream header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_OPEN_FSTREAM_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_OPEN_FSTREAM_HPP

#include <filesystem>
#include <iosfwd>

#include <fmt/core.h>

#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/type_traits.hpp"

namespace cath::common {
	namespace detail {

		/// \brief Function used to implement open_ifstream and open_ofstream overloads that take an fstream
		///
		/// This:
		/// - sets the exceptions to throw on failbit or badbit
		/// - tries to open the file
		/// - throws a runtime_error_exception if the open failed
		///
		/// \param prm_fstream  The fstream to open
		/// \param prm_filename The file to open
		/// \param prm_mode     The mode with which to open the file
		template <typename fstream_t>
		void open_existing_fstream_impl( fstream_t &                      prm_fstream,
		                                 const ::std::filesystem::path &  prm_filename,
		                                 const ::std::ios_base::openmode &prm_mode ) {
			const bool is_reading = is_same_modulo_cvref_v<fstream_t, ::std::ifstream>;
			static_assert( is_reading || is_same_modulo_cvref_v<fstream_t, ::std::ofstream>,
			               "fstream_t must be either ::std::ifstream or ::std::ofstream" );

			prm_fstream.exceptions( ::std::ios::badbit | ( is_reading ? ::std::ios::goodbit : ::std::ios::failbit ) );
			try {
				prm_fstream.open( prm_filename, prm_mode );
			} catch ( const std::exception &ex ) {
				BOOST_THROW_EXCEPTION(
				  cath::common::runtime_error_exception( ::fmt::format( R"(Cannot open file "{}" for {} [{}]: {})",
				                                                        prm_filename.string(),
				                                                        is_reading ? "reading" : "writing",
				                                                        ex.what(),
				                                                        std::strerror( errno ) ) ) );
			};

			if ( !prm_fstream ) {
				BOOST_THROW_EXCEPTION( runtime_error_exception( ::fmt::format( R"(Cannot open file "{}" for {}: {})",
				                                                               prm_filename.string(),
				                                                               is_reading ? "reading" : "writing",
				                                                               strerror( errno ) ) ) );
			}
		}

		/// \brief Function used to implement open_ifstream and open_ofstream overloads that don't take an fstream
		///
		/// This:
		/// - constructs the fstream
		/// - sets the exceptions to throw on failbit or badbit
		/// - tries to open the file
		/// - throws a runtime_error_exception if the open failed
		/// - returns the fstream
		///
		/// \param prm_filename The file to open
		/// \param prm_mode     The mode with which to open the file
		template <typename fstream_t>
		fstream_t open_fstream_impl( const ::std::filesystem::path &prm_filename, const ::std::ios_base::openmode &prm_mode ) {
			fstream_t the_stream;
			open_existing_fstream_impl( the_stream, prm_filename, prm_mode );
			return the_stream;
		}

	} // namespace detail

	void open_ifstream( ::std::ifstream &,
	                    const ::std::filesystem::path &,
	                    const ::std::ios_base::openmode & = ::std::ios_base::in );

	void open_ofstream( ::std::ofstream &,
	                    const ::std::filesystem::path &,
	                    const ::std::ios_base::openmode & = ::std::ios_base::out );

	::std::ifstream open_ifstream( const ::std::filesystem::path &, const ::std::ios_base::openmode & = ::std::ios_base::in );

	::std::ofstream open_ofstream( const ::std::filesystem::path &, const ::std::ios_base::openmode & = ::std::ios_base::out );

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_FILE_OPEN_FSTREAM_HPP
