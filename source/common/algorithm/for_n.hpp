/// \file
/// \brief The for_n() header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_FOR_N_H
#define _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_FOR_N_H

#include <boost/core/ignore_unused.hpp>

#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp17/invoke.hpp"

#include <cstddef>

namespace cath {
	namespace common {

		/// \brief Invoke the specified callable the specified number of times
		template <typename Fn>
		void for_n(const size_t  &arg_n, ///< The number of time to invoke the callable
		           Fn           &&arg_fn ///< The callable to invoke the specified number of times
		           ) {
			for (const auto &x : indices( arg_n ) ) {
				boost::ignore_unused( x );
				invoke( std::forward<Fn>( arg_fn ) );
			}
		}

	} // namespace common
} // namespace cath

#endif
