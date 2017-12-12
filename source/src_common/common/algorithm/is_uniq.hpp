/// \file
/// \brief The is_uniq header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_IS_UNIQ_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_IS_UNIQ_H

#include <boost/range/algorithm/adjacent_find.hpp>
#include <boost/range/concepts.hpp>

#include "common/cpp14/cbegin_cend.hpp"

namespace cath {
	namespace common {

		/// \brief Template function is_uniq()
		///
		/// Returns true iff no to adjacent elements evaluate as equal
		///
		/// Iterator-based algorithm checking whether std::adjacent_find() returns end().
		///
		/// \pre ForwardItr is a model of the Forward Iterator Concept
		template<typename ForwardItr>
		inline bool is_uniq(const ForwardItr &arg_begin, ///< TODOCUMENT
		                    const ForwardItr &arg_end    ///< TODOCUMENT
		                    ) {
			return ( std::adjacent_find( arg_begin, arg_end ) == arg_end );
		}

		/// \brief Template function is_uniq()
		///
		/// Range-based version of the common::is_uniq algorithm, above
		///
		/// \pre ForwardRange is a model of the ForwardRangeConcept
		template <typename ForwardRange>
		inline bool is_uniq(const ForwardRange &arg_rng ///< A single-pass input range
		                    ) {
			BOOST_RANGE_CONCEPT_ASSERT((boost::ForwardRangeConcept<ForwardRange>));
			return is_uniq(
				common::cbegin( arg_rng ),
				common::cend  ( arg_rng )
			);
		}

	} // namespace common
} // namespace cath

#endif
