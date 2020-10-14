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

#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/algorithm/sort_copy.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::numeric_cast;

/// \brief Randomly split the integers [ 0 .. prm_num_instances ) into two vectors, with prm_fraction_in_first (eg 0.5) of them in the first
size_vec_size_vec_pair cath::common::random_split(mt19937       &prm_rng,              ///< TODOCUMENT
                                                  const size_t  &prm_num_instances,    ///< TODOCUMENT
                                                  const double  &prm_fraction_in_first ///< TODOCUMENT
                                                  ) {
	if ( prm_fraction_in_first < 0.0 || prm_fraction_in_first > 1.0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot create random split with fraction in first outside range [0, 1]"));
	}

	auto the_indices = copy_build<size_vec>( indices( prm_num_instances ) );

//	random_shuffle( the_indices, prm_rng );
	shuffle( begin( the_indices ), end( the_indices ), prm_rng );

	const auto num_in_first  = numeric_cast<size_t>( round( prm_fraction_in_first * numeric_cast<double>( prm_num_instances ) ) );
	const auto cut_point_itr = next( cath::common::cbegin( the_indices ), numeric_cast<ptrdiff_t>( num_in_first ) );

	return make_pair(
		sort_copy( size_vec{ cath::common::cbegin( the_indices ), cut_point_itr                 } ),
		sort_copy( size_vec{ cut_point_itr,                   cath::common::cend( the_indices ) } )
	);
}
