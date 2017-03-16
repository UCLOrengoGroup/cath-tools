/// \file
/// \brief The files_equal class header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_TEST_PREDICATE_FILES_EQUAL_H
#define _CATH_TOOLS_SOURCE_COMMON_TEST_PREDICATE_FILES_EQUAL_H

#include <boost/filesystem/path.hpp>
#include <boost/test/test_tools.hpp>

#include "common/test_predicate/istreams_equal.hpp"

namespace cath {

	/// \brief TODOCUMENT
	class files_equal final {
	private:
		/// \brief TODOCUMENT
		bool overwrite_diff_expected_with_got;

		/// \brief TODOCUMENT
		str_size_type diff_half_width;

	public:
		explicit files_equal(const bool &,
		                     const str_size_type & = istreams_equal::DEFAULT_DIFF_HALF_WIDTH);
		explicit files_equal(const str_size_type & = istreams_equal::DEFAULT_DIFF_HALF_WIDTH);

		boost::test_tools::predicate_result operator()(const boost::filesystem::path &,
		                                               const boost::filesystem::path &) const;

		static const std::string FILENAME_NAME_PREFIX;
	};
} // namespace cath

#define BOOST_WARN_FILES_EQUAL(                 S1, S2 )   BOOST_WARN(    ( cath::files_equal(      ) ( ( (S1) ), ( (S2) ) ) ) )
#define BOOST_CHECK_FILES_EQUAL(                S1, S2 )   BOOST_CHECK(   ( cath::files_equal(      ) ( ( (S1) ), ( (S2) ) ) ) )
#define BOOST_REQUIRE_FILES_EQUAL(              S1, S2 )   BOOST_REQUIRE( ( cath::files_equal(      ) ( ( (S1) ), ( (S2) ) ) ) )

#define BOOST_WARN_FILES_EQUAL_OR_OVERWRITE(    S1, S2 )   BOOST_WARN(    ( cath::files_equal( true ) ( ( (S1) ), ( (S2) ) ) ) )
#define BOOST_CHECK_FILES_EQUAL_OR_OVERWRITE(   S1, S2 )   BOOST_CHECK(   ( cath::files_equal( true ) ( ( (S1) ), ( (S2) ) ) ) )
#define BOOST_REQUIRE_FILES_EQUAL_OR_OVERWRITE( S1, S2 )   BOOST_REQUIRE( ( cath::files_equal( true ) ( ( (S1) ), ( (S2) ) ) ) )

#endif
