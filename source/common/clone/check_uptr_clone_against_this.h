/// \file
/// \brief The check_uptr_clone_against_this header

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

#ifndef CHECK_UPTR_CLONE_AGAINST_THIS_H_INCLUDED
#define CHECK_UPTR_CLONE_AGAINST_THIS_H_INCLUDED

#include <boost/core/ignore_unused.hpp>

#include <cassert>
#include <memory>

namespace cath {
	namespace common {

		/// \brief Standard approach to achieving a virtual copy-ctor
		template <typename B, typename D>
		inline std::unique_ptr<B> check_uptr_clone_against_this(std::unique_ptr<B>  arg_clone, ///< TODOCUMENT
		                                                        const D            &arg_clonee ///< TODOCUMENT
		                                                        ) {
			assert( typeid( *arg_clone ) == typeid( arg_clonee ) );
#ifdef NDEBUG
			boost::ignore_unused( arg_clonee );
#endif
			return arg_clone;
		}

	}
}

#endif
