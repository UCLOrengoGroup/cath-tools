/// \file
/// \brief The istreams_equal class definitions

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

#include "istreams_equal.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/test/predicate/detail/strings_equal.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::test;

using ::boost::test_tools::predicate_result;
using ::std::istream;
using ::std::istreambuf_iterator;
using ::std::string;

/// \brief Ctor for istreams_equal
istreams_equal::istreams_equal(const str_size_type &prm_diff_half_width ///< TODOCUMENT
                               ) : diff_half_width(prm_diff_half_width) {
}

/// \brief TODOCUMENT
predicate_result istreams_equal::operator()(istream       &prm_istream1, ///< TODOCUMENT
                                            const string  &prm_name1,    ///< TODOCUMENT
                                            istream       &prm_istream2, ///< TODOCUMENT
                                            const string  &prm_name2     ///< TODOCUMENT
                                            ) const {
	// Suck the two istreams into strings
	return test::detail::strings_equal(
		string{ ( istreambuf_iterator<char>( prm_istream1 ) ), istreambuf_iterator<char>() },
		prm_name1,
		string{ ( istreambuf_iterator<char>( prm_istream2 ) ), istreambuf_iterator<char>() },
		prm_name2,
		diff_half_width
	);
}
