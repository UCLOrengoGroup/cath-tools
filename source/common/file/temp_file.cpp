/// \file
/// \brief The temp_file class definitions

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

#include "temp_file.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>

#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"

#include <iostream> // *** TEMPORARY ***

using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace std;

using boost::none;

/// \brief Constructor for temp_file.
temp_file::temp_file(const string &arg_filename_pattern ///< A pattern for the filename to create with % symbols for characters to be altered (eg ".%%%%-%%%%-%%%%-%%%%.pml")
                     ) : filename( ! arg_filename_pattern.empty() ? opt_path( temp_filename_of_basename_pattern( arg_filename_pattern ) ) : none ) {
}

/// \brief A function to construct a temporary filename from a pattern.
path temp_file::temp_filename_of_basename_pattern(const string &arg_filename_pattern ///< A pattern for the filename to create with % symbols for characters to be altered (eg ".%%%%-%%%%-%%%%-%%%%.pml")
                                                  ) {
	if (path(arg_filename_pattern).has_parent_path()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("temp_file pattern \"" + arg_filename_pattern + "\" should just be a basename, not a full path"));
	}
	return unique_path( temp_directory_path() / arg_filename_pattern );
}

/// \brief Destructor for temp_file
temp_file::~temp_file() noexcept {
	try {
		if ( filename ) {
			if ( exists( *filename ) ) {
				remove( *filename );
			}
		}
	}
	catch (...) {
	}
}
/// \brief TODOCUMENT
const opt_path & temp_file::get_opt_filename() const {
	return filename;
}

/// \brief TODOCUMENT
bool cath::has_filename(const temp_file &arg_temp_file ///< TODOCUMENT
                        ) {
	return static_cast<bool>( arg_temp_file.get_opt_filename() );
}
/// \brief TODOCUMENT
path cath::get_filename(const temp_file &arg_temp_file ///< TODOCUMENT
                        ) {
	return *arg_temp_file.get_opt_filename();
}