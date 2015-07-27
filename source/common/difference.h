/// \file
/// \brief The difference header

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

#ifndef DIFFERENCE_H_INCLUDED
#define DIFFERENCE_H_INCLUDED

#include <algorithm>

namespace cath {
	namespace common {

    	/// \brief Return the difference of two values
    	///
    	/// This is helpful for comparing two values of an unsigned integral type because
    	/// the natural approach of `abs( val_a - val_b )` can hit problems if both are unsigned
    	/// and val_b is bigger than val_a and silly answers can result.
    	///
    	/// That problem doesn't even necessarily get picked up by the compiler. Hence this is a better approach.
    	template <typename T>
    	T difference(const T &arg_value_a, ///< The first value to be compared
    	             const T &arg_value_b  ///< The second value to be compared
    	             ) {
    		return std::max( arg_value_a, arg_value_b ) - std::min( arg_value_a, arg_value_b );
    	}

	}
}

#endif
