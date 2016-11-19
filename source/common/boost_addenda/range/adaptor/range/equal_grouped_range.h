/// \file
/// \brief The equal_grouped_range class header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_RANGE_EQUAL_GROUPED_RANGE_H
#define _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_RANGE_EQUAL_GROUPED_RANGE_H

#include "common/boost_addenda/range/adaptor/iterator/equal_group_itr.h"

namespace cath {
	namespace common {

		/// \brief A range wrapper that dereferences to subranges that are the equivalent groups of elements in the original
		///        range (in the correct order). If the client uses the ctor flag to indicate that this can assume
		///        the original code is sorted in ascending order, then this code needn't do as many comparisons.
		///
		/// For a const equal_grouped_range, use a const RNG type.
		template <typename RNG>
		class equal_grouped_range final : public boost::iterator_range<equal_group_itr<RNG>> {
		private:
			/// \brief The equal_group_itr that does the actual hard work of implementing the equivalent grouping
			using equal_grouped_iterator = equal_group_itr<RNG>;

			/// \brief A convenience type-alias for the iterator_range through which this is implemented
			using super                       = boost::iterator_range<equal_grouped_iterator>;

		public:
			template <typename FN>
			equal_grouped_range(const RNG &,
			                    FN = std::not_equal_to<range_value_t<RNG>>() );
		 };

		/// \brief Ctor from a range and a flag indicating whether it can be assumed that range is sorted
		template <typename RNG>
		template <typename FN>
		equal_grouped_range<RNG>::equal_grouped_range(const RNG &arg_range,           ///< The range over which to apply the equal_grouped_range
		                                              FN         arg_unequal_function ///< TODOCUMENT
		                                              ) : super(
		                                                  	equal_grouped_iterator( std::begin( arg_range ), std::end( arg_range ), arg_unequal_function ),
		                                                  	equal_grouped_iterator( std::end  ( arg_range ), std::end( arg_range ), arg_unequal_function )
		                                                  ) {
		}
	} // namespace common

} // namespace cath

#endif
