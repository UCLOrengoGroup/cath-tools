/// \file
/// \brief The path_step definitions

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
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

#include "path_step.h"

#include <boost/numeric/conversion/cast.hpp>

#include "alignment/alignment.h"
#include "alignment/pair_alignment.h"

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
using namespace std;

using boost::numeric_cast;

/// \brief TODOCUMENT
const path_step_vec path_step_helper::ALL_PATH_STEPS = { path_step::ALIGN_PAIR,
                                                         path_step::INSERT_INTO_FIRST,
                                                         path_step::INSERT_INTO_SECOND };

/// \brief TODOCUMENT
///
/// \relates path_step
ostream & cath::align::detail::operator<<(ostream         &arg_os,       ///< TODOCUMENT
                                          const path_step &arg_path_step ///< TODOCUMENT
                                          ) {
	switch (arg_path_step) {
		case ( path_step::ALIGN_PAIR         ) : { arg_os << "path_step::ALIGN_PAIR"        ; break; }
		case ( path_step::INSERT_INTO_FIRST  ) : { arg_os << "path_step::INSERT_INTO_FIRST" ; break; }
		case ( path_step::INSERT_INTO_SECOND ) : { arg_os << "path_step::INSERT_INTO_SECOND"; break; }
	}
	return arg_os;
}

/// \brief TODOCUMENT
///
/// \relates path_step
size_size_pair cath::align::detail::indices_of_path_step(const path_step &arg_path_step ///< TODOCUMENT
                                                         ) {
	switch(arg_path_step) {
		case( path_step::ALIGN_PAIR ) : {
			return make_pair( 1, 1 );
		}
		case( path_step::INSERT_INTO_FIRST ) : {
			return make_pair( 1, 0 );
		}
		case( path_step::INSERT_INTO_SECOND ) : {
			return make_pair( 0, 1 );
		}
	}
	return make_pair( 0, 0 ); //
}

/// \brief TODOCUMENT
///
/// \relates path_step
size_size_pair cath::align::detail::indices_of_point_after_path_step(const path_step &arg_path_step, ///< TODOCUMENT
                                                                     const size_t    &arg_index_a,   ///< TODOCUMENT
                                                                     const size_t    &arg_index_b    ///< TODOCUMENT
                                                                     ) {
	const size_size_pair path_step_indices = indices_of_path_step( arg_path_step );
	return make_pair(
		arg_index_a + path_step_indices.first,
		arg_index_b + path_step_indices.second
	);
}

/// \brief TODOCUMENT
///
/// \relates path_step
map<path_step, size_size_pair> cath::align::detail::indices_of_point_by_path_step(const size_t &arg_index_a, ///< TODOCUMENT
                                                                                  const size_t &arg_index_b  ///< TODOCUMENT
                                                                                  ) {
	map<path_step, size_size_pair> results;
	for (const path_step &the_path_step : path_step_helper::ALL_PATH_STEPS) {
		results[ the_path_step ] = indices_of_point_after_path_step(
			the_path_step,
			arg_index_a,
			arg_index_b
		);
	}
	return results;
}

/// \brief TODOCUMENT
///
/// \relates path_step
///
/// \relates alignment
void cath::align::detail::append_path_step_to_pair_alignment_from_point(alignment            &arg_alignment, ///< TODOCUMENT
                                                                        const path_step      &arg_path_step, ///< TODOCUMENT
                                                                        const size_size_pair &arg_position   ///< TODOCUMENT
                                                                        ) {
	const auto value_a = arg_position.first;
	const auto value_b = arg_position.second;
	switch ( arg_path_step ) {
		case ( path_step::ALIGN_PAIR         ) : { append_position_both( arg_alignment, value_a, value_b ); break; }
		case ( path_step::INSERT_INTO_FIRST  ) : { append_position_a   ( arg_alignment, value_a          ); break; }
		case ( path_step::INSERT_INTO_SECOND ) : { append_position_b   ( arg_alignment,          value_b ); break; }
	}
	return ;
}
