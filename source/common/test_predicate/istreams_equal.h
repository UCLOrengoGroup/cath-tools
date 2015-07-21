/// \file
/// \brief The istreams_equal class header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Binaries project and then tweaked, eg namespaced in cath)
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

#ifndef ISTREAMS_EQUAL_H_INCLUDED
#define ISTREAMS_EQUAL_H_INCLUDED

#include <boost/test/test_tools.hpp>

#include "common/type_aliases.h"

#include <string>

namespace cath {

	/// \brief TODOCUMENT
	///
	/// Note that the operator() will empty both istreams
	class istreams_equal final {
	private:
		str_size_type diff_half_width;

		static str_size_type index_of_first_difference(const std::string &,
		                                               const std::string &);

	public:
		istreams_equal(const str_size_type &arg_diff_half_width = DEFAULT_DIFF_HALF_WIDTH);

		boost::test_tools::predicate_result operator()(std::istream &,
		                                               const std::string &,
		                                               std::istream &,
		                                               const std::string &) const;

		static const str_size_type DEFAULT_DIFF_HALF_WIDTH;
	};
}

#define BOOST_WARN_ISTREAMS_EQUAL(    I1, S1, I2, S2 )   BOOST_WARN(    ( istreams_equal() ( ( (I1) ), ( (S1) ), ( (I2) ), ( (S2) ) ) ) )
#define BOOST_CHECK_ISTREAMS_EQUAL(   I1, S1, I2, S2 )   BOOST_CHECK(   ( istreams_equal() ( ( (I1) ), ( (S1) ), ( (I2) ), ( (S2) ) ) ) )
#define BOOST_REQUIRE_ISTREAMS_EQUAL( I1, S1, I2, S2 )   BOOST_REQUIRE( ( istreams_equal() ( ( (I1) ), ( (S1) ), ( (I2) ), ( (S2) ) ) ) )

#endif
