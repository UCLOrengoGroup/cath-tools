/// \file
/// \brief The random_split header

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

#ifndef RANDOM_SPLIT_H_INCLUDED
#define RANDOM_SPLIT_H_INCLUDED

#include "common/type_aliases.h"

#include <random>

namespace cath {
	namespace common {

		size_vec_size_vec_pair random_split(std::mt19937 &,
		                                    const size_t &,
		                                    const double &);

	}
}

#endif

