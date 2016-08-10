/// \file
/// \brief The sets_are_disjoint() header

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

#ifndef SETS_ARE_DISJOINT_H_INCLUDED
#define SETS_ARE_DISJOINT_H_INCLUDED

#include <boost/range/algorithm/set_algorithm.hpp>

#include "common/boost_addenda/range/range_concept_type_aliases.h"

#include <vector>

namespace cath {
	namespace common {

		/// \brief Return whether the two set ranges are disjoint (ie have no elements in common)
		///
		/// \tparam RNG1 is a model of the SinglePassRangeConcept
		/// \tparam RNG2 is a model of the SinglePassRangeConcept
		///
		/// \pre RNG1 and RNG2 have the same value type
		///
		/// Note that this is not the most efficient because:
		///  * (a) it could (but does not) stop as soon as it finds any common value
		///  * (b) it goes to the trouble of building vector of the intersection
		///
		/// Still, don't prematurely optimise this until a profiler demonstrates the need because:
		///  (a) it's likely to be fairly quick as is
		///  (b) this is the soft of thing that might end up in other libraries (Boost, std::, range-v3)
		///  (c) rewriting could be quite involved and could involve duplicating quite a bit of std:: code
		///
		/// \todo Consider adding an iterator version
		template <typename RNG1, typename RNG2>
		bool sets_are_disjoint(const RNG1 &arg_range_1, ///< TODOCUMENT
		                       const RNG2 &arg_range_2  ///< TODOCUMENT
		                       ) {
			// Grab the two ranges' value types
			using value_type1 = range_value_t<RNG1>;
			using value_type2 = range_value_t<RNG2>;
			static_assert(
				std::is_same<value_type1, value_type2 >::value,
				"sets_are_disjoint() cannot operator on two ranges with different value types"
			);

			// Build a vector of the values in the intersection
			//
			/// \todo Write a set_intersection_build() and use it here
			std::vector<value_type1> the_intersection;
			boost::range::set_intersection(
				arg_range_1,
				arg_range_2,
				back_inserter( the_intersection )
			);

			// Return whether the intersection is empty
			return the_intersection.empty();
		}
	}
}

#endif
