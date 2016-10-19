/// \file
/// \brief The random_split header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_RANDOM_SPLIT_H
#define _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_RANDOM_SPLIT_H

#include "common/type_aliases.h"

#include <random>

namespace cath {
	namespace common {

		size_vec_size_vec_pair random_split(std::mt19937 &,
		                                    const size_t &,
		                                    const double &);

	} // namespace common
} // namespace cath

#endif

