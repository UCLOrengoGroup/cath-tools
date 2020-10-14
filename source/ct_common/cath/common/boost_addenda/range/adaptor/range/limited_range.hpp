/// \file
/// \brief The limited_range class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_RANGE_LIMITED_RANGE_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_RANGE_LIMITED_RANGE_HPP

#include "cath/common/algorithm/copy_build.hpp" // ***** TEMPORARY *****
#include "cath/common/boost_addenda/range/adaptor/iterator/limit_itr.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"

namespace cath {
	namespace common {

		/// \brief A range wrapper that limits to the first n elements of the underlying range.
		///
		/// For a const limited_range, use a const RNG type.
		///
		/// Invariants:
		///  * The client must preserve the validity of the original range and its one-past-end iterator throughout the limited_range's lifetime
		template <typename RNG>
		class limited_range final : public boost::iterator_range<limit_itr<RNG>> {
		private:
			/// \brief The limit_itr that does the actual hard work of implementing the limit iterator
			using limited_iterator = limit_itr<RNG>;

			/// \brief A convenience type-alias for the iterator_range through which this is implemented
			using super            = boost::iterator_range<limited_iterator>;

		public:
			limited_range(const RNG &,
			              const size_t &);
		 };

		/// \brief Ctor from a range
		template <typename RNG>
		limited_range<RNG>::limited_range(const RNG    &prm_range,           ///< The range over which to apply this limited_range
		                                  const size_t &prm_max_num_elements ///< The maximum number of elements
		                                  ) : super(
		                                      	limited_iterator( cath::common::cbegin( prm_range ), cath::common::cend  ( prm_range ), prm_max_num_elements ),
		                                      	limited_iterator( cath::common::cend  ( prm_range ), cath::common::cend  ( prm_range ), prm_max_num_elements )
		                                      ) {
		}

	} // namespace common
} // namespace cath

#endif
