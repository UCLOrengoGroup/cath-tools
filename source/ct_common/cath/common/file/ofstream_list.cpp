/// \file
/// \brief The ofstream_list class definitions

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

#include "ofstream_list.hpp"

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/file/open_fstream.hpp"

#include <fstream>

using namespace ::cath;
using namespace ::cath::common;

using ::boost::filesystem::path;
using ::std::ofstream;
using ::std::ostream;

/// \brief Ctor from a special ostream and a flag to be used to indicate when output should be sent to that ostream
ofstream_list::ofstream_list(ostream    &prm_standard_outstream,     ///< A special ostream (often stdout) to which output can be sent
                             const path &prm_standard_outstream_flag ///< A flag to be used to indicate when output should be sent to the special ostream
                             ) : standard_outstream     { prm_standard_outstream      },
                                 standard_outstream_flag{ prm_standard_outstream_flag } {
}

/// \brief Open the specified list of paths and add them to the outputs
ostream_ref_vec ofstream_list::open_ofstreams(const path_vec &prm_paths ///< The paths to open
                                              ) {
	return transform_build<ostream_ref_vec>(
		prm_paths,
		[&] (const path &the_path) -> ostream_ref {
			if ( the_path == standard_outstream_flag && standard_outstream ) {
				return *standard_outstream;
			}
			else {
				ofstreams.emplace_back();
				open_ofstream( ofstreams.back(), the_path );
				return { ofstreams.back() };
			}
		}
	);
}

/// \brief Get the flag, which is used to indicate when output should be sent to the special ostream
const path & ofstream_list::get_flag() const {
	return standard_outstream_flag;
}

/// \brief Close all ofstreams
void ofstream_list::close_all() {
	for (ofstream &the_ofstream: ofstreams) {
		the_ofstream.close();
	}
}

/// \brief Open a single path and add it to the outputs in the specified ofstream_list
///
/// \relates ofstream_list
ostream_ref cath::common::open_ofstream(ofstream_list &prm_ofstreams, ///< The ofstream_list in which to open the path
                                        const path    &prm_path       ///< The path to open
                                        ) {
	ostream_ref_vec ostream_refs = prm_ofstreams.open_ofstreams( { prm_path } );
	if ( ostream_refs.size() != 1 ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Did not get back one ostream from opening one file"));
	}
	ostream_ref result = std::move( ostream_refs.back() );
	return result;
}
