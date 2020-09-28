/// \file
/// \brief The append header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_APPEND_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_APPEND_HPP

#include <iterator>
#include <utility>

#include <boost/range/algorithm/copy.hpp>

namespace cath {
	namespace common {

		/// \brief Append the specified range to the specified container
		///
		/// Note: previously used insert from boost/range/algorithm_ext/insert.hpp but
		/// that appeared to sometimes cause problems detected by Clang Asan.
		///
		/// \param prm_container The container to which the range should be appended
		/// \param prm_suffix    The range to append to the container
		template <typename Cont, typename Rng>
		inline Cont &append( Cont &prm_container, Rng &&prm_suffix ) {
			::boost::range::copy( ::std::forward<Rng>( prm_suffix ), ::std::back_inserter( prm_container ) );
			return prm_container;
		}

		/// \brief Append the specified range to a copy of the specified container
		///
		/// \param prm_container The container to which the range should be appended
		/// \param prm_suffix    The range to append to the container
		template <typename Cont, typename Rng>
		inline Cont append_copy( Cont prm_container, Rng &&prm_suffix ) {
			append( prm_container, ::std::forward<Rng>( prm_suffix ) );
			return prm_container;
		}

	} // namespace common
} // namespace cath

#endif
