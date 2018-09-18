/// \file
/// \brief The coord_list class definitions

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

#include "coord_list.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/numeric.hpp>

#include "common/cpp14/cbegin_cend.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "structure/geometry/coord.hpp"

#include <iostream> // ***** TEMPORARY *****
#include <numeric>

using namespace cath::common;
using namespace cath::geom;
using namespace std;

using boost::accumulate;
using boost::numeric_cast;
using boost::inner_product;

/// \brief Constructor for coord_list from a vector of coords.
coord_list::coord_list(coord_vec prm_coords ///< The vector of coords with which this coord_list should be initialised
                       ) : coords{ std::move( prm_coords ) } {
}

/// \brief Standard reserve() method for reserving memory for more coords
void coord_list::reserve(const size_t &prm_size ///< The number of coords that this coord list should have adequate memory to store
                         ) {
	coords.reserve( prm_size );
}

/// \brief Standard empty() to return whether the coord list is empty (ie contains zero coords)
bool coord_list::empty() const noexcept {
	return coords.empty();
}

/// \brief Standard size() operator to return the number of coords in the coord list
size_t coord_list::size() const {
	return coords.size();
}

/// \brief Standard push_back method to add a coord the back of a list
void coord_list::push_back(const coord &prm_coord ///< The coord to push on the back of the list
                           ) {
	return coords.push_back( prm_coord );
}

/// \brief Non-const-overload of standard subscript operator
coord & coord_list::operator[](const size_t &prm_index ///< The index of the coord in the coord_list to access
                               ) {
	return coords[ prm_index ];
}

/// \brief Const-overload of standard subscript operator
const coord & coord_list::operator[](const size_t &prm_index ///< The index of the coord in the coord_list to access
                                     ) const {
	return coords[ prm_index ];
}

/// \brief Subtract a coord add to all entries in the coord list
void coord_list::operator+=(const coord &prm_coord ///< The coord to add to all entries in the coord list
                            ) {
	for (coord &coord_loop : coords) {
		coord_loop += prm_coord;
	}
}

/// \brief Subtract a coord from all entries in the coord list
void coord_list::operator-=(const coord &prm_coord ///< The coord to subtract from all entries in the coord list
                            ) {
	for (coord &coord_loop : coords) {
		coord_loop -= prm_coord;
	}
}

/// \brief Standard non-const begin() operator to provide range access
coord_list::iterator coord_list::begin() {
	return std::begin( coords );
}
/// \brief Standard non-const end() operator to provide range access
coord_list::iterator coord_list::end() {
	return std::end( coords );
}
/// \brief Standard const begin() operator to provide range access
coord_list::const_iterator coord_list::begin() const {
	return common::cbegin( coords );
}
/// \brief Standard const end() operator to provide range access
coord_list::const_iterator coord_list::end() const {
	return common::cend( coords );
}

/// \brief Flatten a vector of one coord_list
///
/// \todo Make a `flattened` range adaptor and then see whether all calls to this can be replaced with that.
///
/// \relates coord_list
coord_list cath::geom::flatten_coord_lists(const coord_list_vec &prm_coord_lists ///< TODOCUMENT
                                           ) {
	coord_list result;
	for (const coord_list &the_coord_list : prm_coord_lists) {
		for (const coord &the_coord : the_coord_list) {
			result.push_back( the_coord );
		}
	}
	return result;
}

/// \brief Calculate the sum coord of a list of coords
///
/// \relates coord_list
coord cath::geom::sum(const coord_list &prm_coord_list ///< The list of coords to sum
                      ) {
	if ( prm_coord_list.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot calculate centre of gravity for empty coord_list"));
	}
	return accumulate( prm_coord_list, coord::ORIGIN_COORD );
}

/// \brief Grab the number of coords in each coord_list and throw if they're zero or not equal
///
/// \relates coord_list
size_t cath::geom::check_non_empty_and_equal_size(const coord_list &prm_coord_list_1, ///< The first  list of coords to compare
                                                  const coord_list &prm_coord_list_2  ///< The second list of coords to compare
                                                  ) {
	const auto num_coords_1 = prm_coord_list_1.size();
	const auto num_coords_2 = prm_coord_list_2.size();
	if (num_coords_1 != num_coords_2) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("coord_lists must be of equal size for calc_rmsd()"));
	}
	if (num_coords_1 == 0) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("coord_lists must be non-empty for calc_rmsd()"));
	}
	return num_coords_1;
}

/// \brief Calculate the centre of geometry from a list of coords
///
/// \relates coord_list
coord cath::geom::centre_of_gravity(const coord_list &prm_coord_list ///< The list of coords for which the centre_of_gravity should be calculated
                                    ) {
	return sum( prm_coord_list ) / numeric_cast<double>( prm_coord_list.size() );
}

/// \brief Calculate the mean deviation between two lists of coords
///
/// \relates coord_list
double cath::geom::calc_mean_deviation(const coord_list &prm_coord_list_1, ///< The first  list of coords to compare
                                       const coord_list &prm_coord_list_2  ///< The second list of coords to compare
                                       ) {
	// \brief Grab the number of coords in each coord_list and throw if they're zero or not equal
	const auto num_coords = check_non_empty_and_equal_size(
		prm_coord_list_1,
		prm_coord_list_2
	);

	// Calculate the sum of the squared distances between pairs of corresponding coords
	const auto total_deviation = inner_product(
		prm_coord_list_1,
		prm_coord_list_2,
		0.0,
		plus<double>(),
//		[] (const coord &x, const coord &y) { cerr << "Distance between " << x << " and " << y << " is " << distance_between_points( x, y ) << endl; return distance_between_points( x, y ); }
		[] (const coord &x, const coord &y) { return distance_between_points( x, y ); }
	);

	// Return the total deviation by the number of coords to get the MD
	return total_deviation / numeric_cast<double>( num_coords );
}

/// \brief Calculate the RMSD between two lists of coords
///
/// \relates coord_list
double cath::geom::calc_rmsd(const coord_list &prm_coord_list_1, ///< The first  list of coords to compare
                             const coord_list &prm_coord_list_2  ///< The second list of coords to compare
                             ) {
	// \brief Grab the number of coords in each coord_list and throw if they're zero or not equal
	const auto num_coords = check_non_empty_and_equal_size(
		prm_coord_list_1,
		prm_coord_list_2
	);

	// Calculate the sum of the squared distances between pairs of corresponding coords
	const auto total_squared_deviation = inner_product(
		prm_coord_list_1,
		prm_coord_list_2,
		0.0,
		plus<double>(),
		[] (const coord &x, const coord &y) { return squared_distance_between_points( x, y ); }
	);

	// Divide the SD by the number of coords to get the MSD
	const auto mean_squared_deviation = total_squared_deviation / numeric_cast<double>( num_coords );

	// Return the square root to give the RMSD
	return sqrt( mean_squared_deviation );
}

/// \brief Basic insertion operator to output a rough summary of an coord_list to an ostream
///
/// \relates coord_list
ostream & cath::geom::operator<<(ostream          &prm_os,        ///< The ostream to which to output the summary
                                 const coord_list &prm_coord_list ///< The coord_list to summarise
                                 ) {
	prm_os << "coord_list[\n";
	for (const coord &curr_coord : prm_coord_list) {
		prm_os << "\t" << curr_coord << "\n";
	}
	prm_os << "]";
	return prm_os;
}
