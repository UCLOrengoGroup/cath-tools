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

#ifndef CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_PREDICATE_ISTREAM_AND_FILE_EQUAL_HPP
#define CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_PREDICATE_ISTREAM_AND_FILE_EQUAL_HPP

#include <filesystem>
#include <string>

#include <boost/test/test_tools.hpp>

#include "cath/test/predicate/bootstrap_mode.hpp"
#include "cath/test/predicate/istreams_equal.hpp"

namespace cath::test {

	/// \brief TODOCUMENT
	///
	/// Note that the operator() will empty the istream
	class istream_and_file_equal final {
	  private:
		/// \brief When to bootstrap the test file (ie replace it with the "got" content if mismatching)
		bootstrap_mode bootstrapping = DEFAULT_BOOTSTRAPPING;

		/// \brief TODOCUMENT
		str_size_type diff_half_width = istreams_equal::DEFAULT_DIFF_HALF_WIDTH;

		/// \brief Default value for when to bootstrap the test file (ie replace it with the "got" content if mismatching)
		static constexpr bootstrap_mode DEFAULT_BOOTSTRAPPING = bootstrap_mode::IF_ENV;

	  public:
		explicit istream_and_file_equal( const bootstrap_mode &,
		                                 const str_size_type & = istreams_equal::DEFAULT_DIFF_HALF_WIDTH );
		explicit istream_and_file_equal( const str_size_type & = istreams_equal::DEFAULT_DIFF_HALF_WIDTH );
		boost::test_tools::predicate_result operator()( std::istream &, const std::string &, const ::std::filesystem::path & ) const;
	};

} // namespace cath::test

// clang-format off
#define BOOST_WARN_ISTREAM_AND_FILE_EQUAL(    I1, S1, F2 ) BOOST_WARN(    ( cath::test::istream_and_file_equal( ) ( ( (I1) ), ( (S1) ), ( (F2) ) ) ) )
#define BOOST_CHECK_ISTREAM_AND_FILE_EQUAL(   I1, S1, F2 ) BOOST_CHECK(   ( cath::test::istream_and_file_equal( ) ( ( (I1) ), ( (S1) ), ( (F2) ) ) ) )
#define BOOST_REQUIRE_ISTREAM_AND_FILE_EQUAL( I1, S1, F2 ) BOOST_REQUIRE( ( cath::test::istream_and_file_equal( ) ( ( (I1) ), ( (S1) ), ( (F2) ) ) ) )
// clang-format on

#endif // CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_PREDICATE_ISTREAM_AND_FILE_EQUAL_HPP
