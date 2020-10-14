/// \file
/// \brief The istreams_equal class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_TEST_TEST_PREDICATE_ISTREAMS_EQUAL_HPP
#define _CATH_TOOLS_SOURCE_SRC_TEST_TEST_PREDICATE_ISTREAMS_EQUAL_HPP

#include <boost/test/test_tools.hpp>

#include "cath/common/type_aliases.hpp"

#include <string>

namespace cath {
	namespace test {

		/// \brief TODOCUMENT
		///
		/// Note that the operator() will empty both istreams
		class istreams_equal final {
		private:
			/// \brief TODOCUMENT
			str_size_type diff_half_width;

			static str_size_type index_of_first_difference(const std::string &,
			                                               const std::string &);

		public:
			explicit istreams_equal(const str_size_type & = DEFAULT_DIFF_HALF_WIDTH);

			boost::test_tools::predicate_result operator()(std::istream &,
			                                               const std::string &,
			                                               std::istream &,
			                                               const std::string &) const;

			/// \brief The default half-width used when displaying any differences
			static constexpr str_size_type DEFAULT_DIFF_HALF_WIDTH = 50;
		};

	} // namespace test
} // namespace cath

#define BOOST_WARN_ISTREAMS_EQUAL(    I1, S1, I2, S2 )   BOOST_WARN(    ( cath::test::istreams_equal() ( ( (I1) ), ( (S1) ), ( (I2) ), ( (S2) ) ) ) )
#define BOOST_CHECK_ISTREAMS_EQUAL(   I1, S1, I2, S2 )   BOOST_CHECK(   ( cath::test::istreams_equal() ( ( (I1) ), ( (S1) ), ( (I2) ), ( (S2) ) ) ) )
#define BOOST_REQUIRE_ISTREAMS_EQUAL( I1, S1, I2, S2 )   BOOST_REQUIRE( ( cath::test::istreams_equal() ( ( (I1) ), ( (S1) ), ( (I2) ), ( (S2) ) ) ) )

#endif
