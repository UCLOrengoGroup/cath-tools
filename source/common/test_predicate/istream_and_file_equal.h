/// \file
/// \brief The istream_and_file_equal class header

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

#ifndef ISTREAM_AND_FILE_EQUAL_H_INCLUDED
#define ISTREAM_AND_FILE_EQUAL_H_INCLUDED

#include <boost/filesystem/path.hpp>
#include <boost/test/test_tools.hpp>

#include "common/test_predicate/istreams_equal.h"

#include <string>

namespace cath {

	/// \brief TODOCUMENT
	///
	/// Note that the operator() will empty the istream
	class istream_and_file_equal final {
	private:
		/// \brief TODOCUMENT
		bool overwrite_diff_expected_with_got;

		/// \brief TODOCUMENT
		str_size_type diff_half_width;

	public:
		istream_and_file_equal(const bool &,
		                       const str_size_type &arg_diff_half_width = istreams_equal::DEFAULT_DIFF_HALF_WIDTH);
		explicit istream_and_file_equal(const str_size_type &arg_diff_half_width = istreams_equal::DEFAULT_DIFF_HALF_WIDTH);
		boost::test_tools::predicate_result operator()(std::istream &,
		                                               const std::string &,
		                                               const boost::filesystem::path &) const;
	};
}

#define BOOST_WARN_ISTREAM_AND_FILE_EQUAL(                 I1, S1, F2 )   BOOST_WARN(    ( istream_and_file_equal(      ) ( ( (I1) ), ( (S1) ), ( (F2) ) ) ) )
#define BOOST_CHECK_ISTREAM_AND_FILE_EQUAL(                I1, S1, F2 )   BOOST_CHECK(   ( istream_and_file_equal(      ) ( ( (I1) ), ( (S1) ), ( (F2) ) ) ) )
#define BOOST_REQUIRE_ISTREAM_AND_FILE_EQUAL(              I1, S1, F2 )   BOOST_REQUIRE( ( istream_and_file_equal(      ) ( ( (I1) ), ( (S1) ), ( (F2) ) ) ) )

#define BOOST_WARN_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE(    I1, S1, F2 )   BOOST_WARN(    ( istream_and_file_equal( true ) ( ( (I1) ), ( (S1) ), ( (F2) ) ) ) )
#define BOOST_CHECK_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE(   I1, S1, F2 )   BOOST_CHECK(   ( istream_and_file_equal( true ) ( ( (I1) ), ( (S1) ), ( (F2) ) ) ) )
#define BOOST_REQUIRE_ISTREAM_AND_FILE_EQUAL_OR_OVERWRITE( I1, S1, F2 )   BOOST_REQUIRE( ( istream_and_file_equal( true ) ( ( (I1) ), ( (S1) ), ( (F2) ) ) ) )

#endif
