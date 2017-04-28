/// \file
/// \brief The remove_itrs_from_range() header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_COPY_BUILD_H
#define _CATH_TOOLS_SOURCE_COMMON_ALGORITHM_COPY_BUILD_H

#include <boost/range/irange.hpp>

#include "common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/size_t_literal.hpp"

namespace cath {
	namespace common {

		/// \brief Remove the elements at the specified positions from the specified range
		///
		/// \returns The new end of the range (a la std::range)
		template <typename Rng,
		          typename RmvRng>
		range_iterator_t<Rng> remove_itrs_from_range(Rng          &arg_range,  ///< The range from which the elements should be removed
		                                             const RmvRng &arg_removes ///< A range of (possibly const_iterator) iterators into arg_range, which are to be removed
		                                             ) {
			if ( arg_removes.empty() ) {
				return std::end( arg_range );
			}

			auto write_itr = next(
				std::begin( arg_range ),
				distance( common::cbegin( arg_range ), arg_removes.front() )
			);

			const auto move_range_fn = [&] (auto begin_itr, const auto end_itr) {
				while ( begin_itr != end_itr ) {
					*write_itr = std::move( *begin_itr );
					++write_itr;
					++begin_itr;
				}
			};
			for (const size_t remove_ctr : boost::irange( 1_z, arg_removes.size() ) ) {
				move_range_fn(
					next( arg_removes[ remove_ctr - 1 ] ),
					      arg_removes[ remove_ctr     ]
				);
			}

			move_range_fn(
				next( arg_removes.back() ),
				common::cend( arg_range )
			);

			return write_itr;
		}

	} // namespace common
} // namespace cath
#endif
