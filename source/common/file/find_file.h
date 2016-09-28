/// \file
/// \brief The find_file class header

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

#ifndef FIND_FILE_H_INCLUDED
#define FIND_FILE_H_INCLUDED

#include <boost/filesystem/path.hpp>

#include "common/path_type_aliases.h"

#include <string>

namespace cath {
	namespace common {

		boost::filesystem::path find_file(const path_vec &,
		                                  const std::string &);

	}
}

#endif
