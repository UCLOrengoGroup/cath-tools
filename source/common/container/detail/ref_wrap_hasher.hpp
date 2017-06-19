/// \file
/// \brief The ref_wrap_hasher header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_CONTAINER_DETAIL_REF_WRAP_HASHER_H
#define _CATH_TOOLS_SOURCE_COMMON_CONTAINER_DETAIL_REF_WRAP_HASHER_H

#include "common/cpp14/cbegin_cend.hpp"

#include <functional>
#include <type_traits>

namespace cath { namespace common { namespace detail { template <typename T> class ref_wrap_uom_wrap; } } }

namespace cath {
	namespace common {

		namespace detail {

			/// \brief Hasher for a ref_wrap_uom_wrap<T>
			template <typename T>
			struct ref_wrap_hasher final {
				/// \brief The function operator that performs the hash on the T value
				///        using std::hash<decay_t<T>>
				size_t operator()(const ref_wrap_uom_wrap<T> &arg_value
				                  ) const {
					return std::hash<std::decay_t<T>>{}( arg_value.get() );
				}
			};
		}

	} // namespace common
} // namespace cath

#endif
