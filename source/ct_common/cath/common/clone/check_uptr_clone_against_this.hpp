/// \file
/// \brief The check_uptr_clone_against_this header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CLONE_CHECK_UPTR_CLONE_AGAINST_THIS_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CLONE_CHECK_UPTR_CLONE_AGAINST_THIS_HPP

#include <cassert>
#include <memory>

namespace cath::common {

	/// \brief Standard approach to achieving a virtual copy-ctor
	///
	/// \param prm_clone  TODOCUMENT
	/// \param prm_clonee TODOCUMENT
	template <typename B, typename D>
	inline std::unique_ptr<B> check_uptr_clone_against_this( std::unique_ptr<B> prm_clone, [[maybe_unused]] const D &prm_clonee ) {
		assert( typeid( *prm_clone ) == typeid( prm_clonee ) );
		return prm_clone;
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CLONE_CHECK_UPTR_CLONE_AGAINST_THIS_HPP
