/// \file
/// \brief The sorted_indices header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_SORTED_INDICES_H
#define _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_SORTED_INDICES_H

#include <boost/range/algorithm/sort.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/boost_addenda/range/indices.hpp"

namespace cath {
	namespace common {

			/// \brief Build a container of indices [0, prm_size), sorted by the specified less-than function
			///
			/// \pre prm_size > 0
			template <typename Cont = std::vector<size_t>, typename Fn>
			Cont sorted_indices(const size_t  &prm_size,    ///< The number of consecutive indices (starting from 0) to sort
			                    Fn           &&prm_less_fn  ///< The less-than function with which to sort the indices
			                    ) {
				auto sorting_indices = common::copy_build<Cont>( common::indices( prm_size ) );
				boost::range::sort(
					sorting_indices,
					std::forward<Fn>( prm_less_fn )
				);
				return sorting_indices;
			}

	} // namespace common
} // namespace cath

#endif
