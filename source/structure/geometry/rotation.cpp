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

#include "rotation.h"

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "exception/invalid_argument_exception.h"
#include "exception/runtime_error_exception.h"
#include "structure/geometry/angle.h"
#include "structure/geometry/coord.h"
#include "structure/geometry/coord_list.h"

#include <iostream> // **** TEMPORARY ****

using namespace boost::math;
using namespace boost::test_tools;
using namespace cath::common;
using namespace cath::geom;
using namespace std;

using boost::lexical_cast;

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

const double & rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS() {
	static const double default_tolerance_for_rotation_closeness_checks(0.000001);
	return default_tolerance_for_rotation_closeness_checks;
}

/// \brief Non-member equality operator for the rotation class
///
/// \relates rotation
bool cath::geom::operator==(const rotation &arg_rot_a, ///< TODOCUMENT
                            const rotation &arg_rot_b  ///< TODOCUMENT
                            ) {
	for (size_t new_row_ctr = 0; new_row_ctr < coord::NUM_DIMS; ++new_row_ctr) {
		for (size_t new_col_ctr = 0; new_col_ctr < coord::NUM_DIMS; ++new_col_ctr) {
			if ( difference( arg_rot_a.get_value( new_row_ctr, new_col_ctr ), arg_rot_b.get_value( new_row_ctr, new_col_ctr ) ) > rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS() ) {
				return false;
			}
		}
	}

	return true;
}

/// \brief TODOCUMENT
///
/// \relates rotation
std::ostream & cath::geom::operator<<(std::ostream   &arg_os,      ///< TODOCUMENT
                                      const rotation &arg_rotation ///< TODOCUMENT
                                      ) {
	arg_os << "rotation[";
	for (size_t new_row_ctr = 0; new_row_ctr < coord::NUM_DIMS; ++new_row_ctr) {
		arg_os << ((new_row_ctr == 0) ? "" : "; ");
		for (size_t new_col_ctr = 0; new_col_ctr < coord::NUM_DIMS; ++new_col_ctr) {
			arg_os << ((new_col_ctr == 0) ? "" : ", ");
			arg_os << right << setw( 7 ) << arg_rotation.get_value(new_row_ctr, new_col_ctr);
		}
	}
	arg_os << "]";
	return arg_os;
}

/// \brief Return whether two rotations are very close to each other.
///
/// \relates rotation
///
/// This compares each pair of values and uses rotation::TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS to determine "close".
bool cath::geom::are_close(const rotation &arg_rotn_1, ///< TODOCUMENT
                           const rotation &arg_rotn_2  ///< TODOCUMENT
                           ) {
	for (size_t new_row_ctr = 0; new_row_ctr < coord::NUM_DIMS; ++new_row_ctr) {
		for (size_t new_col_ctr = 0; new_col_ctr < coord::NUM_DIMS; ++new_col_ctr) {
			if ( difference( arg_rotn_1.get_value( new_row_ctr, new_col_ctr ), arg_rotn_2.get_value( new_row_ctr, new_col_ctr ) ) > rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS() ) {
				return false;
			}
		}
	}
	return true;
}

/// \brief Finds the rotation that rotates the first coord on the x-axis and the second on the x-y plane, with positive y
///
/// \relates rotation
rotation cath::geom::rotation_to_x_axis_and_x_y_plane(const coord &arg_a, ///< TODOCUMENT
                                                      const coord &arg_b  ///< TODOCUMENT
                                                      ) {
	assert(not_zero(arg_a));
	assert(not_zero(arg_b));
	const coord normalised_a          = normalise_copy(arg_a);
	assert(not_zero(cross_product(normalised_a, arg_b)));
	const coord orthogonal_to_a_and_b = normalise_copy(cross_product(normalised_a, arg_b));
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
                                   const double &arg_tolerance ///< The maximum difference that should be tolerated between the input values and the result before an exception is thrown
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
//	for (size_t col_ctr = 0; col_ctr < coord::NUM_DIMS; ++col_ctr) {
//		for (size_t row_ctr = 0; row_ctr < coord::NUM_DIMS; ++row_ctr) {
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
//			<< ( difference( values[ value_index_ctr_temp ], new_rotation.get_value( row_ctr, col_ctr ) ) > arg_tolerance )
//			<< endl;
//			++value_index_ctr_temp;
//		}
//	}

	size_t value_index_ctr = 0;
	for (size_t row_ctr = 0; row_ctr < coord::NUM_DIMS; ++row_ctr) {
		for (size_t col_ctr = 0; col_ctr < coord::NUM_DIMS; ++col_ctr) {
			if ( difference( values[ value_index_ctr ], new_rotation.get_value( row_ctr, col_ctr ) ) > arg_tolerance ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception(
					"Unable to create rotation that is within tolerance "
					+ boost::lexical_cast<std::string>( arg_tolerance )
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
rotation cath::geom::tidy_copy(const rotation &arg_rotation, ///< The rotation matrix to tidy up
                               const double   &arg_tolerance ///< The maximum difference that should be tolerated between the input values and the result before an exception is thrown
                               ) {
	return tidy_rotation(
		arg_rotation.get_value<0, 0>(), arg_rotation.get_value<0, 1>(), arg_rotation.get_value<0, 2>(),
		arg_rotation.get_value<1, 0>(), arg_rotation.get_value<1, 1>(), arg_rotation.get_value<1, 2>(),
		arg_rotation.get_value<2, 0>(), arg_rotation.get_value<2, 1>(), arg_rotation.get_value<2, 2>(),
		arg_tolerance
	);
}



