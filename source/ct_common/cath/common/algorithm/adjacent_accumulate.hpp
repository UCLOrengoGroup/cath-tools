/// \file
/// \brief The adjacent_accumulate() header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_ADJACENT_ACCUMULATE_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_ADJACENT_ACCUMULATE_HPP

#include <boost/range/concepts.hpp>

#include "cath/common/cpp14/cbegin_cend.hpp"

namespace cath {
	namespace common {

		/// \brief Accumulate the results of executing a binary function on the pairs of adjacent values in a range
		///
		/// \pre Container        is a model of the Mutable_Container concept
		/// \pre SinglePassRange1 is a model of the SinglePassRangeConcept
		template <typename SinglePassItr,
		          typename Value,
		          typename BinaryOperation>
		inline Value adjacent_accumulate(SinglePassItr        prm_begin, ///< Iterator to the beginning of a single-pass input range
		                                 const SinglePassItr &prm_end,   ///< Iterator to the end of a single-pass input range
		                                 Value                prm_init,  ///< An initial value
		                                 BinaryOperation      prm_fun    ///< A binary function to execute on the adjacent pairs of elements of the range
		                                 ) {
			while ( std::next( prm_begin ) != prm_end ) {
				prm_init += prm_fun( *prm_begin, *std::next( prm_begin ) );
				++prm_begin;
			}
			return prm_init;
		}

		/// \brief Accumulate the results of executing a binary function on the pairs of adjacent values in a range
		///
		/// WARNING: adjacent_difference() is a bit rubbish - it copies the first element of the input range into
		///          the output range before doing calculating the adjacent_differences. In practice, this means
		///          the output range must have the same value type and the input range.
		///
		/// \sa adjacent_accumulate() generate_n_build() random_sample_n_build() sort_adjacent_difference() sort_uniq_adjacent_difference() transform_build() uniq_adjacent_difference()
		///
		/// \pre Container        is a model of the Mutable_Container concept
		/// \pre SinglePassRange1 is a model of the SinglePassRangeConcept
		template <typename SinglePassRange,
		          typename Value,
		          typename BinaryOperation>
		inline Value adjacent_accumulate(const SinglePassRange &prm_rng, ///< A single-pass input range
		                                 Value                  prm_init, ///< An initial value
		                                 BinaryOperation        fun       ///< A binary function to execute on the adjacent pairs of elements of the range
		                                 ) {
			// Check that SinglePassRange is a SinglePassRangeConcept
			BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const SinglePassRange > ));

			return adjacent_accumulate(
				common::cbegin( prm_rng ),
				common::cend  ( prm_rng ),
				std::move( prm_init ),
				fun
			);
		}


	} // namespace common
} // namespace cath
#endif
