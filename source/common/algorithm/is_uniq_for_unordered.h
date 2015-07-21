/// \file
/// \brief The is_uniq_for_unordered header

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

#ifndef IS_UNIQ_FOR_UNORDERED_H_INCLUDED
#define IS_UNIQ_FOR_UNORDERED_H_INCLUDED

#include <boost/range.hpp>

namespace cath {
	namespace common {
		/// \brief TODOCUMENT
		template <typename ForwardItr>
		inline bool is_uniq_for_unordered(const ForwardItr &arg_begin, ///< TODOCUMENT
		                                  const ForwardItr &arg_end    ///< TODOCUMENT
		                                  ) {
			for (ForwardItr itr( arg_begin ); itr != arg_end; ++itr ) {
				if ( count( itr, arg_end, *itr ) > 1 ) {
					return false;
				}
			}
			return true;
		}

		/// \brief TODOCUMENT
		template <typename ForwardRange>
		inline bool is_uniq_for_unordered(const ForwardRange &arg_rng ///< TODOCUMENT
		                                  ) {
			BOOST_RANGE_CONCEPT_ASSERT((boost::ForwardRangeConcept<ForwardRange>));
			return is_uniq_for_unordered(
				common::cbegin( arg_rng ),
				common::cend  ( arg_rng )
			);
		}
	}
}

#endif
