/// \file
/// \brief The path_or_istream class definitions

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

#include "path_or_istream.hpp"

using namespace cath::common;

using boost::filesystem::path;
using boost::none;
using std::istream;

/// \brief Ctor from a special istream and a flag to be used to indicate when input should be read from that istream
path_or_istream::path_or_istream(istream    &arg_istream,               ///< A special istream (often stdin) from which input can be read
                                 const path &arg_standard_instream_flag ///< A flag to be used to indicate when input should be read from the special istream
                                 ) : standard_instream      { arg_istream                },
                                     standard_instream_flag { arg_standard_instream_flag } {
}

/// \brief Open the specified list of paths and add them to the inputs
path_or_istream & path_or_istream::set_path(const path &arg_file ///< The path to open
                                            ) {
	if ( input_file_stream ) {
		input_file_stream->close();
	}
	if ( arg_file != get_flag() || ! standard_instream ) {
		input_file_stream.emplace();
		open_ifstream( *input_file_stream, arg_file );
	}
	return *this;
}

/// \brief Get the flag, which is used to indicate when input should be read from the special istream
const path & path_or_istream::get_flag() const {
	return standard_instream_flag;
}

/// \brief Close the ifstream, if any
path_or_istream & path_or_istream::close() {
	if ( input_file_stream ) {
		input_file_stream->close();
		input_file_stream = none;
	}
	return *this;
}

/// \brief Get the active istream
istream & path_or_istream::get_istream() {
	return input_file_stream ? static_cast<istream &>( *input_file_stream )
	                         : standard_instream->get();
}

/// \brief Return whether the specified file will actually be used by the specified
///        path_or_istream (ie is not the special flag) but is missing
///
/// \relates path_or_istream
bool cath::common::file_is_missing(const path_or_istream &arg_path_or_istream, ///< The path_or_istream containing the special flag to indicate whether the input should be read from an istream rather than a file
                                   const path            &arg_file             ///< The file in question
                                   ) {
	return (
		arg_file != arg_path_or_istream.get_flag()
		&&
		! exists( arg_file )
	);
}
