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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_DETAIL_REF_WRAP_HASHER_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_DETAIL_REF_WRAP_HASHER_HPP

#include <functional>

#include "cath/common/type_traits.hpp"

// clang-format off
namespace cath::common::detail { template <typename T> class ref_wrap_uom_wrap; }
// clang-format on

namespace cath::common::detail {

	/// \brief Hasher for a ref_wrap_uom_wrap<T>
	template <typename T>
	struct ref_wrap_hasher final {
		/// \brief The function operator that performs the hash on the T value
		///        using std::hash<remove_cvref_t<T>>
		size_t operator()(const ref_wrap_uom_wrap<T> &prm_value
		                  ) const {
			return std::hash<common::remove_cvref_t<T>>{}( prm_value.get() );
		}
	};

} // namespace cath::common::detail

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CONTAINER_DETAIL_REF_WRAP_HASHER_HPP
