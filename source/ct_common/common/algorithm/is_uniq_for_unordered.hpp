/// \file
/// \brief The is_uniq_for_unordered header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_IS_UNIQ_FOR_UNORDERED_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_IS_UNIQ_FOR_UNORDERED_HPP

#include <boost/range.hpp>

#include "common/cpp14/cbegin_cend.hpp"

#include <algorithm>

namespace cath {
	namespace common {
		/// \brief TODOCUMENT
		template <typename ForwardItr>
		inline bool is_uniq_for_unordered(const ForwardItr &prm_begin, ///< TODOCUMENT
		                                  const ForwardItr &prm_end    ///< TODOCUMENT
		                                  ) {
			for (ForwardItr itr( prm_begin ); itr != prm_end; ++itr ) {
				if ( std::count( itr, prm_end, *itr ) > 1 ) {
					return false;
				}
			}
			return true;
		}

		/// \brief TODOCUMENT
		template <typename ForwardRange>
		inline bool is_uniq_for_unordered(const ForwardRange &prm_rng ///< TODOCUMENT
		                                  ) {
			BOOST_RANGE_CONCEPT_ASSERT((boost::ForwardRangeConcept<ForwardRange>));
			return is_uniq_for_unordered(
				common::cbegin( prm_rng ),
				common::cend  ( prm_rng )
			);
		}
	} // namespace common
} // namespace cath

#endif
