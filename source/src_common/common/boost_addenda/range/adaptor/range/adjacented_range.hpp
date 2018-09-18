/// \file
/// \brief The adjacented_range class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_RANGE_ADJACENTED_RANGE_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_RANGE_ADJACENTED_RANGE_HPP

#include "common/boost_addenda/range/adaptor/iterator/adjacent_itr.hpp"
#include "common/cpp14/cbegin_cend.hpp"

namespace cath {
	namespace common {

		/// \brief A range wrapper that dereferences to the consecutive pairs of adjacent elements.
		///
		/// For a const adjacented_range, use a const RNG type.
		///
		/// \todo If there is need, this could possibly be extended to allow compile-time specification of the number
		///       of elements to return at once (via a tuple)
		///
		/// Invariants:
		///  * The client must preserve the validity of the original range and its one-past-end iterator throughout the adjacented_range's lifetime
		template <typename RNG>
		class adjacented_range final : public boost::iterator_range<adjacent_itr<RNG>> {
		private:
			/// \brief The adjacent_itr that does the actual hard work of implementing the adjacent iterator
			using adjacented_iterator = adjacent_itr<RNG>;

			/// \brief A convenience type-alias for the iterator_range through which this is implemented
			using super               = boost::iterator_range<adjacented_iterator>;

		public:
			explicit adjacented_range(const RNG &);
		 };

		/// \brief Ctor from a range
		template <typename RNG>
		adjacented_range<RNG>::adjacented_range(const RNG &prm_range ///< The range over which to apply this adjacented_range
		                                        ) : super(
		                                            	adjacented_iterator( cath::common::cbegin( prm_range ), cath::common::cend  ( prm_range ) ),
		                                            	adjacented_iterator( cath::common::cend  ( prm_range ), cath::common::cend  ( prm_range ) )
		                                            ) {
		}

	} // namespace common
} // namespace cath

#endif
