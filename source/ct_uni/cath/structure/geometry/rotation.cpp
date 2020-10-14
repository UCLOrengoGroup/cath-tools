/// \file
/// \brief The rotation class definitions

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

#include "rotation.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/structure/geometry/angle.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/coord_list.hpp"

#include <iostream> // **** TEMPORARY ****

using namespace ::cath::common;
using namespace ::cath::geom;
using namespace ::std;

using ::boost::algorithm::any_of;
using ::boost::lexical_cast;
using ::boost::property_tree::ptree;

/// \brief The identity rotation
const rotation & rotation::IDENTITY_ROTATION() {
	static const rotation identity_rotation(1.0, 0.0, 0.0,
	                                        0.0, 1.0, 0.0,
	                                        0.0, 0.0, 1.0);
	return identity_rotation;
}

/// \brief TODOCUMENT
const rotation & rotation::ROTATE_X_TO_Y_TO_Z_TO_X() {
	static const rotation rotate_x_to_y_to_z_to_x(0.0, 0.0, 1.0,
	                                              1.0, 0.0, 0.0,
	                                              0.0, 1.0, 0.0);
	return rotate_x_to_y_to_z_to_x;
}

/// \brief TODOCUMENT
const rotation & rotation::ROTATE_X_TO_Z_TO_Y_TO_X() {
	static const rotation rotate_x_to_z_to_y_to_x(0.0, 1.0, 0.0,
	                                              0.0, 0.0, 1.0,
	                                              1.0, 0.0, 0.0);
	return rotate_x_to_z_to_y_to_x;
}

constexpr double rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS;

/// \brief Non-member equality operator for the rotation class
///
/// \relates rotation
bool cath::geom::operator==(const rotation &prm_rot_a, ///< TODOCUMENT
                            const rotation &prm_rot_b  ///< TODOCUMENT
                            ) {
	for (const size_t &new_row_ctr : indices( coord::NUM_DIMS ) ) {
		for (const size_t &new_col_ctr : indices( coord::NUM_DIMS ) ) {
			if ( difference( prm_rot_a.get_value( new_row_ctr, new_col_ctr ), prm_rot_b.get_value( new_row_ctr, new_col_ctr ) ) > rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS ) {
				return false;
			}
		}
	}

	return true;
}

/// \brief TODOCUMENT
///
/// \relates rotation
std::ostream & cath::geom::operator<<(std::ostream   &prm_os,      ///< TODOCUMENT
                                      const rotation &prm_rotation ///< TODOCUMENT
                                      ) {
	prm_os << "rotation[";
	for (const size_t &new_row_ctr : indices( coord::NUM_DIMS ) ) {
		prm_os << ((new_row_ctr == 0) ? "" : "; ");
		for (const size_t &new_col_ctr : indices( coord::NUM_DIMS ) ) {
			prm_os << ((new_col_ctr == 0) ? "" : ", ");
			prm_os << right << setw( 7 ) << prm_rotation.get_value(new_row_ctr, new_col_ctr);
		}
	}
	prm_os << "]";
	return prm_os;
}

/// \brief Return whether two rotations are very close to each other.
///
/// \relates rotation
///
/// This compares each pair of values and uses rotation::TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS to determine "close".
bool cath::geom::are_close(const rotation &prm_rotn_1, ///< TODOCUMENT
                           const rotation &prm_rotn_2  ///< TODOCUMENT
                           ) {
	for (const size_t &new_row_ctr : indices( coord::NUM_DIMS ) ) {
		for (const size_t &new_col_ctr : indices( coord::NUM_DIMS ) ) {
			if ( difference( prm_rotn_1.get_value( new_row_ctr, new_col_ctr ), prm_rotn_2.get_value( new_row_ctr, new_col_ctr ) ) > rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS ) {
				return false;
			}
		}
	}
	return true;
}

/// \brief Finds the rotation that rotates the first coord on the x-axis and the second on the x-y plane, with positive y
///
/// \relates rotation
rotation cath::geom::rotation_to_x_axis_and_x_y_plane(const coord &prm_a, ///< TODOCUMENT
                                                      const coord &prm_b  ///< TODOCUMENT
                                                      ) {
	assert(not_zero(prm_a));
	assert(not_zero(prm_b));
	const coord normalised_a          = normalise_copy(prm_a);
	assert(not_zero(cross_product(normalised_a, prm_b)));
	const coord orthogonal_to_a_and_b = normalise_copy(cross_product(normalised_a, prm_b));
	const coord b_perpendicular       = normalise_copy(cross_product(orthogonal_to_a_and_b, normalised_a));

	return rotation(
			     normalised_a.get_x(),          normalised_a.get_y(),          normalised_a.get_z(),
		      b_perpendicular.get_x(),       b_perpendicular.get_y(),       b_perpendicular.get_z(),
		orthogonal_to_a_and_b.get_x(), orthogonal_to_a_and_b.get_y(), orthogonal_to_a_and_b.get_z()
	);
}

/// \brief Tidy up the values for a rotation matrix if they are not quite precise enough
///
/// \relates rotation
///
/// Note that the arguments use row-major order so that they can be laid out like the matrix.
rotation cath::geom::tidy_rotation(const double &val_00, const double &val_01, const double &val_02,
                                   const double &val_10, const double &val_11, const double &val_12,
                                   const double &val_20, const double &val_21, const double &val_22,
                                   const double &prm_tolerance ///< The maximum difference that should be tolerated between the input values and the result before an exception is thrown
                                   ) {
	/// Attempt to build a new matrix like the original but of higher precision
	const rotation new_rotation = rotation_to_x_axis_and_x_y_plane(
		coord( val_00, val_01, val_02 ),
		coord( val_10, val_11, val_12 )
	);

	// Put the original values in a vector
	const doub_vec values = { val_00, val_01, val_02,
	                          val_10, val_11, val_12,
	                          val_20, val_21, val_22 };

	/// Check that each of the values in the new matrix is adequately close to the equivalent input value

//	cerr << endl;
//	size_t value_index_ctr_temp = 0;
//	for (const size_t &col_ctr : indices( coord::NUM_DIMS ) ) {
//		for (const size_t &row_ctr : indices( coord::NUM_DIMS ) ) {
//			cerr << "Row "
//			<< lexical_cast<string>( row_ctr )
//			<< " and column "
//			<< lexical_cast<string>( col_ctr )
//			<< "; specified "
//			<< lexical_cast<string>( values[value_index_ctr_temp] )
//			<< "; got "
//			<< lexical_cast<string>( new_rotation.get_value(row_ctr, col_ctr) )
//			<< ")\t"
//			<< boolalpha
//			<< ( difference( values[ value_index_ctr_temp ], new_rotation.get_value( row_ctr, col_ctr ) ) > prm_tolerance )
//			<< endl;
//			++value_index_ctr_temp;
//		}
//	}

	size_t value_index_ctr = 0;
	for (const size_t &row_ctr : indices( coord::NUM_DIMS ) ) {
		for (const size_t &col_ctr : indices( coord::NUM_DIMS ) ) {
			if ( difference( values[ value_index_ctr ], new_rotation.get_value( row_ctr, col_ctr ) ) > prm_tolerance ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception(
					"Unable to create rotation that is within tolerance "
					+ boost::lexical_cast<std::string>( prm_tolerance )
					+ " to specified elements (failure at row "
					+ boost::lexical_cast<std::string>( row_ctr )
					+ " and column "
					+ boost::lexical_cast<std::string>( col_ctr )
					+ "; specified "
					+ boost::lexical_cast<std::string>( values[value_index_ctr] )
					+ " but got "
					+ boost::lexical_cast<std::string>( new_rotation.get_value(row_ctr, col_ctr) )
					+ ")"
				));
			}
			++value_index_ctr;
		}
	}

	return new_rotation;
}

/// \brief Tidy up the values for a rotation matrix if they are not quite precise enough
///
/// \relates rotation
rotation cath::geom::tidy_copy(const rotation &prm_rotation, ///< The rotation matrix to tidy up
                               const double   &prm_tolerance ///< The maximum difference that should be tolerated between the input values and the result before an exception is thrown
                               ) {
	return tidy_rotation(
		prm_rotation.get_value<0, 0>(), prm_rotation.get_value<0, 1>(), prm_rotation.get_value<0, 2>(),
		prm_rotation.get_value<1, 0>(), prm_rotation.get_value<1, 1>(), prm_rotation.get_value<1, 2>(),
		prm_rotation.get_value<2, 0>(), prm_rotation.get_value<2, 1>(), prm_rotation.get_value<2, 2>(),
		prm_tolerance
	);
}

/// \brief Make a new Boost Property Tree ptree containing a single anonymous entry for the specified rotation's
///        value at the specified row and column
///
/// \relates rotation
ptree cath::geom::detail::make_ptree_of_row_and_col(const rotation &prm_rotation, ///< The rotation containing the entry to be saved
                                                    const size_t   &prm_row,      ///< The row index of the entry to be saved within the specified rotation
                                                    const size_t   &prm_column    ///< The column index of the entry to be saved within the specified rotation
                                                    ) {
	ptree new_ptree;
	new_ptree.put( "", prm_rotation.get_value( prm_row, prm_column ) );
	return new_ptree;
}

/// \brief Make a new Boost Property Tree ptree containing an anonymous array of anonymous entries for the specified rotation's
///        values in the specified row
///
/// \relates rotation
ptree cath::geom::detail::make_ptree_of_row(const rotation &prm_rotation, ///< The rotation containing the row to be saved
                                            const size_t   &prm_row       ///< The row index of the specified rotation to be saved
                                            ) {
	ptree new_ptree;
	new_ptree.push_back( make_pair( "", make_ptree_of_row_and_col( prm_rotation, prm_row, 0 ) ) );
	new_ptree.push_back( make_pair( "", make_ptree_of_row_and_col( prm_rotation, prm_row, 1 ) ) );
	new_ptree.push_back( make_pair( "", make_ptree_of_row_and_col( prm_rotation, prm_row, 2 ) ) );
	return new_ptree;
}

/// \brief Build a rotation from a rotation-populated ptree
///
/// \relates rotation
rotation cath::geom::rotation_from_ptree(const ptree &prm_ptree ///< The ptree from which the rotation should be read
                                         ) {
	if ( prm_ptree.size() != coord::NUM_DIMS || prm_ptree.count( "" ) != coord::NUM_DIMS ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse rotation from property_tree that doesn't have three anonymous entries for the three rows"));
	}
	if ( any_of( prm_ptree, [&] (const pair<const string, ptree> &x) { return ( x.second.size() != coord::NUM_DIMS || x.second.count( "" ) != coord::NUM_DIMS ); } ) ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse rotation from property_tree with a row that doesn't contain three anonymous entries"));
	}
	doub_vec values;
	values.reserve( coord::NUM_DIMS * coord::NUM_DIMS );
	for (const auto &row : prm_ptree) {
		for (const auto &row_col : row.second ) {
			values.push_back( row_col.second.get<double>( "" ) );
		}
	}
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return rotation{ values };
}

/// \brief Save the specified rotation to the specified Boost Property Tree ptree
///
/// \relates rotation
void cath::geom::save_to_ptree(ptree          &prm_ptree,   ///< The ptree to which the rotation should be saved
                               const rotation &prm_rotation ///< The rotation to be saved
                               ) {
	prm_ptree.push_back( make_pair( "", detail::make_ptree_of_row( prm_rotation, 0 ) ) );
	prm_ptree.push_back( make_pair( "", detail::make_ptree_of_row( prm_rotation, 1 ) ) );
	prm_ptree.push_back( make_pair( "", detail::make_ptree_of_row( prm_rotation, 2 ) ) );
}
