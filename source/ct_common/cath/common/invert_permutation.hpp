/// \file
/// \brief The invert_permutation() header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_INVERT_PERMUTATION_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_INVERT_PERMUTATION_H

#include <boost/range/combine.hpp>
#include <boost/range/size.hpp>
#include <boost/tuple/tuple.hpp>

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace common {

		/// \brief Return an inverted copy of the specified permutation
		///
		/// A permutation is a reordered list of the integers from 0 to some number (inclusive)
		///
		/// Example: { 4, 1, 0, 3, 2 } inverts to { 2, 1, 4, 3, 0 } and vice versa
		///
		/// This is useful for converting (an ordering of indices) into
		/// (a lookup from index to rank in the ordering)
		template <typename Cont = std::vector<size_t>, typename Rng>
		Cont invert_permutation(const Rng &prm_range ///< The range containing the permutation to invert
		                        ) {
			const size_t range_size = boost::size( prm_range );
			Cont index_scores( range_size, range_size );
			for (const auto &val : combine( common::indices( range_size ), prm_range ) ) {
				index_scores[ static_cast<size_t>( boost::get<1>( val ) ) ] = boost::get<0>( val );
			}
			return index_scores;
		}

	} // namespace common
} // namespace cath

#endif
