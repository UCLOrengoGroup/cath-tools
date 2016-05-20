/// \file
/// \brief The make_clone header

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

#ifndef MAKE_CLONE_H_INCLUDED
#define MAKE_CLONE_H_INCLUDED

#include <memory>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Template for generating clones from a value
			///
			/// This allows specialisation for types that need to use a method name other than "clone",
			/// if clone() has already been used for producing a clone at an earlier level in the hierarchy.
			template <typename T>
			std::unique_ptr<T> make_clone(const T &arg_value ///< TODOCUMENT
			                              ) {
				return arg_value.clone();
			}
		}
	}
}

#endif
