/// \file
/// \brief The windowed_matrix class definitions

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

#include "windowed_matrix.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/not_implemented_exception.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "common/temp_check_offset_1.hpp"

#include <algorithm>

using namespace cath;
using namespace cath::common;
using namespace std;

/// \brief Private method for converting two indices to a internal index
size_t windowed_matrix::get_internal_index(const size_type &/*prm_index_a*/,
                                           const size_type &/*prm_index_b*/
                                           ) const {
	BOOST_THROW_EXCEPTION(not_implemented_exception("Not implemented matrix indexing yet"));
	return 0;
}

/// \brief Ctor for windowed_matrix
windowed_matrix::windowed_matrix(const size_type &prm_length_a,        ///< The length of the first  entry (or the number of rows    in     matrix)
                                 const size_type &prm_length_b,        ///< The length of the second entry (or the number of columns in the matrix)
                                 const size_type &prm_requested_window ///< The requested width (or equivalently, height) of the window
                                 ) : length_a             ( prm_length_a                ),
                                     length_b             ( prm_length_b                ),
                                     requested_window_size( prm_requested_window        ),
                                     data                 ( prm_length_a * prm_length_b ) {
	BOOST_THROW_EXCEPTION(not_implemented_exception("Not implemented matrix size"));
	check_lengths_and_window_size_are_valid( length_a, length_b, requested_window_size );
}

/// \brief TODOCUMENT
windowed_matrix::size_type windowed_matrix::get_length_a() const {
	return length_a;
}

/// \brief TODOCUMENT
windowed_matrix::size_type windowed_matrix::get_length_b() const {
	return length_b;
}

/// \brief TODOCUMENT
windowed_matrix::size_type windowed_matrix::get_window_size() const {
	return requested_window_size;
}

/// \brief TODOCUMENT
void windowed_matrix::set_value(const size_type &prm_index_a, ///< TODOCUMENT
                                const size_type &prm_index_b, ///< TODOCUMENT
                                const double    &prm_value    ///< TODOCUMENT
                                ) {
	const size_t internal_index = get_internal_index( prm_index_a, prm_index_b );
	data[ internal_index ] = prm_value;
}

/// \brief TODOCUMENT
double windowed_matrix::get_value(const size_type &prm_index_a, ///< TODOCUMENT
                                  const size_type &prm_index_b  ///< TODOCUMENT
                                  ) const {
	const size_t internal_index = get_internal_index( prm_index_a, prm_index_b );
	return data[ internal_index ];
}

