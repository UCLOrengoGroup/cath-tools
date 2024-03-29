/// \file
/// \brief The make_uptr_clone class header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CLONE_MAKE_UPTR_CLONE_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CLONE_MAKE_UPTR_CLONE_HPP

#include <memory>

namespace cath::common {

	/// \brief Make a unique_ptr clone of the specified object using its copy-ctor
	///
	/// This allows do_clone() methods to be implemented in one line:
	///
	/// ~~~~~.cpp
	/// return { make_uptr_clone( *this ) };
	/// ~~~~~
	///
	/// \param prm_clonee The object to be cloned (via its copy-ctor) into a unique_ptr
	template <typename T>
	inline auto make_uptr_clone( const T &prm_clonee ) {
		return ::std::make_unique<T>( prm_clonee );
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CLONE_MAKE_UPTR_CLONE_HPP
