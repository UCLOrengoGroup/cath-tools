/// \file
/// \brief The is_cloneable concept template class

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
///
/// This class is designed for checking that a type X provides a const clone() method
/// whose return type can be used to construct an unique_ptr<X>.
///
/// This allows for: T is a const type but T's clone() returns an unique_ptr to non-const T.
///
/// \todo Test and document whether it also allows for a clone that returns a raw pointer.
///
/// This uses the Boost Concept tools.
/// See the "Creating Concept Checking Classes" page for the Boost Concept Check Library pages
/// (currently http://www.boost.org/libs/concept_check/creating_concepts.htm)

#ifndef IS_CLONEABLE_H_INCLUDED
#define IS_CLONEABLE_H_INCLUDED

#include <boost/concept_check.hpp>

#include "common/clone/detail/make_clone.h"

namespace cath {
	namespace common {
		namespace detail {

			/// \brief TODOCUMENT
			template <typename X>
			struct is_cloneable final {
				BOOST_CONCEPT_USAGE(is_cloneable) {
					std::unique_ptr<X> new_ptr( make_clone( value ) );
				}
				const X &value;
			};

		}
	}
}

#endif
