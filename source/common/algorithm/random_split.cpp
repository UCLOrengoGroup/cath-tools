/// \file
/// \brief The random_split definitions

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

#include "random_split.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/sort_copy.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/size_t_literal.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::irange;
using boost::numeric_cast;
using boost::range::random_shuffle;

/// \brief Randomly split the integers [ 0 .. arg_num_instances ) into two vectors, with arg_fraction_in_first (eg 0.5) of them in the first
size_vec_size_vec_pair cath::common::random_split(mt19937       &arg_rng,              ///< TODOCUMENT
                                                  const size_t  &arg_num_instances,    ///< TODOCUMENT
                                                  const double  &arg_fraction_in_first ///< TODOCUMENT
                                                  ) {
	if ( arg_fraction_in_first < 0.0 || arg_fraction_in_first > 1.0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot create random split with fraction in first outside range [0, 1]"));
	}

	auto indices = copy_build<size_vec>( irange( 0_z, arg_num_instances ) );

//	random_shuffle( indices, arg_rng );
	shuffle( begin( indices ), end( indices ), arg_rng );

	const auto num_in_first  = numeric_cast<size_t>( round( arg_fraction_in_first * numeric_cast<double>( arg_num_instances ) ) );
	const auto cut_point_itr = next( cath::common::cbegin( indices ), numeric_cast<ptrdiff_t>( num_in_first ) );

	return make_pair(
		sort_copy( size_vec{ cath::common::cbegin( indices ), cut_point_itr                 } ),
		sort_copy( size_vec{ cut_point_itr,                   cath::common::cend( indices ) } )
	);
}
