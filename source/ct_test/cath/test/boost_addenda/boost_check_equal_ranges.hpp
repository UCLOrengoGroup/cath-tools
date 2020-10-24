/// \file
/// \brief The BOOST_CHECK_EQUAL_RANGES()header

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

#ifndef _CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_BOOST_ADDENDA_BOOST_CHECK_EQUAL_RANGES_HPP
#define _CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_BOOST_ADDENDA_BOOST_CHECK_EQUAL_RANGES_HPP

#include <boost/range/concepts.hpp>
#include <boost/test/test_tools.hpp>

#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"

namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		///
		/// \todo Consider is it worth making this into a macro so that the test output
		///       provides the line numbers and variable names from the caller's context
		///       rather than from here?
		template <typename RNG1, typename RNG2>
		void BOOST_CHECK_EQUAL_RANGES(const RNG1 &prm_rng_1, ///< TODOCUMENT
		                              const RNG2 &prm_rng_2  ///< TODOCUMENT
		                              ) {
			BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const RNG1 > ));
			BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept< const RNG2 > ));

			BOOST_CHECK_EQUAL_COLLECTIONS(
				common::cbegin( prm_rng_1 ),
				common::cend  ( prm_rng_1 ),
				common::cbegin( prm_rng_2 ),
				common::cend  ( prm_rng_2 )
			);
		}


	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_BOOST_ADDENDA_BOOST_CHECK_EQUAL_RANGES_HPP
