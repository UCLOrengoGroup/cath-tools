/// \file
/// \brief The string_matches_file class header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_TEST_PREDICATE_STRING_MATCHES_FILE_H
#define _CATH_TOOLS_SOURCE_COMMON_TEST_PREDICATE_STRING_MATCHES_FILE_H

#include <boost/filesystem/path.hpp>
#include <boost/test/test_tools.hpp>

#include "common/test_predicate/bootstrap_mode.hpp"
#include "common/test_predicate/istreams_equal.hpp"

#include <string>

namespace cath {
	namespace test {

		/// \brief TODOCUMENT
		///
		/// Note that the operator() will empty the string
		class string_matches_file final {
		private:
			/// \brief When to bootstrap the test file (ie replace it with the "got" content if mismatching)
			bootstrap_mode bootstrapping   = DEFAULT_BOOTSTRAPPING;

			/// \brief TODOCUMENT
			str_size_type  diff_half_width = istreams_equal::DEFAULT_DIFF_HALF_WIDTH;


			/// \brief Default value for when to bootstrap the test file (ie replace it with the "got" content if mismatching)
			static constexpr bootstrap_mode DEFAULT_BOOTSTRAPPING = bootstrap_mode::IF_ENV;

		public:
			/// \brief TODOCUMENT
			static constexpr str_size_type  DEFAULT_DIFF_HALF_WIDTH                  = istreams_equal::DEFAULT_DIFF_HALF_WIDTH;

			string_matches_file() noexcept = default;
			explicit string_matches_file(const bootstrap_mode &,
			                             const str_size_type & = istreams_equal::DEFAULT_DIFF_HALF_WIDTH);
			explicit string_matches_file(const str_size_type &);
			boost::test_tools::predicate_result operator()(const std::string &,
			                                               const boost::filesystem::path &) const;
		};

	} // namespace test
} // namespace cath

#define BOOST_WARN_STRING_MATCHES_FILE(    F1, S1 ) BOOST_WARN   ( ( cath::test::string_matches_file( ) ( ( (F1) ), ( (S1) ) ) ) )
#define BOOST_CHECK_STRING_MATCHES_FILE(   F1, S1 ) BOOST_CHECK  ( ( cath::test::string_matches_file( ) ( ( (F1) ), ( (S1) ) ) ) )
#define BOOST_REQUIRE_STRING_MATCHES_FILE( F1, S1 ) BOOST_REQUIRE( ( cath::test::string_matches_file( ) ( ( (F1) ), ( (S1) ) ) ) )

#endif
