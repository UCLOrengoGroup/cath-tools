/// \file
/// \brief The invert_permutation class definitions

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

#include "invert_permutation.h"

#include <boost/lexical_cast.hpp>

#include "exception/invalid_argument_exception.h"

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;

/// \brief Return an inverted copy of the specified permutation
///
/// A permutation is a reordered list of the integers from 0 to some number (inclusive)
///
/// Example: { 4, 1, 0, 3, 2 } inverts to { 2, 1, 4, 3, 0 } and vice versa
size_vec cath::invert_permutation(const size_vec &arg_permutation ///< The permutation to invert
                                  ) {
	const size_t permutation_size        = arg_permutation.size();
	const size_t not_yet_populated_value = permutation_size + 1;
	size_vec results(permutation_size, not_yet_populated_value);
	for (size_t permutation_from = 0; permutation_from < permutation_size; ++permutation_from) {
		const size_t &permutation_to = arg_permutation[ permutation_from ];
		if ( permutation_to >= permutation_size ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Cannot invert permutation because it is invalid: entry "
				+ lexical_cast<string>( permutation_from )
				+ " has value "
				+ lexical_cast<string>( permutation_to   )
				+ ", which is too large for a permutation of size "
				+ lexical_cast<string>( permutation_size )
			));
		}
		if ( results[permutation_to] != not_yet_populated_value ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Cannot invert permutation because it is invalid: value "
				+ lexical_cast<string>( permutation_to )
				+ " appears multiple times"
			));
		}
		results[permutation_to] = permutation_from;
	}
	return results;
}
