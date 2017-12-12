/// \file
/// \brief The rotation class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_GEOMETRY_ROTATION_H
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_GEOMETRY_ROTATION_H

#include <boost/algorithm/clamp.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/operators.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

#include "common/boost_addenda/range/indices.hpp"
#include "common/difference.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "common/exception/runtime_error_exception.hpp"
#include "common/property_tree/read_from_ptree.hpp"
#include "common/type_aliases.hpp"
#include "structure/geometry/angle.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/coord_list.hpp"

#include <iostream> // ***** TEMPORARY *****
#include <vector>

namespace cath {
	namespace geom {

		/// \brief Represent a three-dimensional rotation matrix
		///
		/// ATM, this class is read-only after construction because this keeps things simple.
		///
		/// The class implements basic post-construction checks that this is a valid rotation by checking that:
		///  - multiplying on either side by the transpose gives a matrix very close to the identity
		///  - the determinant is very close to +1
		///
		/// Since these checks are implemented after every construction, everything else can reasonably assume that a rotation is valid.
		///
		/// \todo Change this to be implemented via some standard matrix/rotation class. Boost?
		///
		/// If you're not very familiar with rotation matrices, it might be worth skimming some reference
		/// (eg http://en.wikipedia.org/wiki/Rotation_matrix)
		/// but here are a few pointers:
		///  - Every rotation in 3D is a rotation of a certain angle around a certain axis
		///  - Rotations are applied to vectors through multiplication, so the vector v rotated by the rotation A is A * v
		///  - A rotation's transpose acts as its inverse (ie it does the opposite rotation back to the original orientation)
		///  - Multiplying rotations A*B gives the rotation found by rotating by B and then by A
		///  - Rotation multiplication is associative (ie it is generally true that A*(B*C) = (A*B)*C
		///  - Rotation multiplication isn't commutative (ie it is not generally true that A*B = B*A)
		class rotation final : private boost::multipliable<rotation>,
		                       private boost::equality_comparable<rotation> {
		private:
			/// \brief TODOCUMENT
			double value_0_0;
			/// \brief TODOCUMENT
			double value_0_1;
			/// \brief TODOCUMENT
			double value_0_2;
			/// \brief TODOCUMENT
			double value_1_0;
			/// \brief TODOCUMENT
			double value_1_1;
			/// \brief TODOCUMENT
			double value_1_2;
			/// \brief TODOCUMENT
			double value_2_0;
			/// \brief TODOCUMENT
			double value_2_1;
			/// \brief TODOCUMENT
			double value_2_2;

	//		/// \brief The rotation matrix
	//		doub_vec_vec matrix;

			/// \brief TODOCUMENT
			double tolerance;

			static void check_index(const size_t &);
			static void check_value(const double &);
			rotation & init_from_vector_with_rotation_checks(const doub_vec &);
			rotation & init_from_vector(const doub_vec &);
			rotation & init_from_values_with_rotation_checks(const double &, const double &, const double &,
			                                                 const double &, const double &, const double &,
			                                                 const double &, const double &, const double &);
			rotation & init_from_values(const double &, const double &, const double &,
			                            const double &, const double &, const double &,
			                            const double &, const double &, const double &);
			rotation & set_value(const size_t &, const size_t &, const double &);
			void check_is_valid_rotation() const;
			double determinant() const;



			template <size_t index>
			void check_index() const;

			template <size_t row_index, size_t col_index>
			void set_value(const double &);

		public:
			rotation(const double &, const double &, const double &,
			         const double &, const double &, const double &,
			         const double &, const double &, const double &,
			         const double & = DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS);
			explicit rotation(const doub_vec &,
			                  const double & = DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS);

			const double & get_value(const size_t &, const size_t &) const;

			template <size_t row_index, size_t col_index>
			const double & get_value() const;

			void operator*=(const rotation &);

			static const rotation & IDENTITY_ROTATION();
			static const rotation & ROTATE_X_TO_Y_TO_Z_TO_X();
			static const rotation & ROTATE_X_TO_Z_TO_Y_TO_X();

			static constexpr double DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS = 0.000001;
		};


		bool operator==(const rotation &, const rotation &);
		std::ostream & operator<<(std::ostream &, const rotation &);

		template<size_t row_index>
		coord get_row(const rotation &);
		template<size_t col_index>
		coord get_column(const rotation &);

		coord get_row(const rotation &,
		              const size_t &);
		coord get_column(const rotation &,
		                 const size_t &);

		bool are_close(const rotation &, const rotation &);
		rotation transpose_copy(const rotation &);
		rotation rotation_between_rotations(const rotation &,
		                                    const rotation &);
		double trace(const rotation &);
		doub_angle angle_of_rotation(const rotation &);
		doub_angle angle_between_rotations(const rotation &,
		                                   const rotation &);

		rotation rotation_to_x_axis_and_x_y_plane(const coord &, const coord &);
		rotation tidy_rotation(const double &, const double &, const double &,
		                       const double &, const double &, const double &,
		                       const double &, const double &, const double &,
		                       const double & = rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS);
		rotation tidy_copy(const rotation &,
		                   const double & = rotation::DEFAULT_TOLERANCE_FOR_ROTATION_CLOSENESS_CHECKS);

		void rotate(const rotation &, coord_list &);
		void rotate(const rotation &, coord &);
		coord_list rotate_copy(const rotation &, const coord_list &);
		coord      rotate_copy(const rotation &, const coord &);

		namespace detail {
			boost::property_tree::ptree make_ptree_of_row_and_col(const rotation &,
			                                                      const size_t &,
			                                                      const size_t &);

			boost::property_tree::ptree make_ptree_of_row(const rotation &,
			                                              const size_t &);
		} // namespace detail

		rotation rotation_from_ptree(const boost::property_tree::ptree &);

		void save_to_ptree(boost::property_tree::ptree &,
		                   const rotation &);

		/// \brief Generate some rotation that has the specified angle of rotation
		///
		/// \relates rotation
		template <typename T>
		rotation rotation_of_angle(const angle<T> &arg_angle /// The angle through which the created rotation should rotate
		                           ) {
			const double angle_in_rads = angle_in_radians( arg_angle );
			return rotation_to_x_axis_and_x_y_plane(
				coord( cos( angle_in_rads ), 0.0, sin( angle_in_rads ) ),
				coord( 0.0, 1.0, 0.0 )
			);
		}

		/// \brief Check that a dimension index is a sensible value
		///        (or else throw an exception)
		inline void rotation::check_index(const size_t &arg_index ///< TODOCUMENT
		                                  ) {
	#ifndef NDEBUG
			if (arg_index >= coord::NUM_DIMS) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Index is out of range for rotation matrix"));
			}
	#else
			boost::ignore_unused( arg_index );
	#endif
		}

		/// \brief Check that a dimension index is a sensible value
		///        (or else throw an exception)
		inline void rotation::check_value(const double &arg_value
		                                  ) {
	#ifndef NDEBUG
			if ( ! boost::math::isfinite( arg_value ) ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Value is not a sensible finite number"));
			}
	#else
			boost::ignore_unused( arg_value );
	#endif
		}

		/// \brief Private initialisation from row-major order vector to be used by multiple ctors
		///        and then basic checks that this is a valid rotation
		inline rotation & rotation::init_from_vector_with_rotation_checks(const doub_vec &arg_vector ///< TODOCUMENT
		                                                                  ) {
			init_from_vector( arg_vector );
			check_is_valid_rotation();
			return *this;
		}

		/// \brief Private initialisation from row-major order vector to be used by multiple ctors
		inline rotation & rotation::init_from_vector(const doub_vec &arg_vector ///< TODOCUMENT
		                                             ) {
			// Sanity check the inputs
			if ( arg_vector.size() != coord::NUM_DIMS * coord::NUM_DIMS ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Rotation matrix cannot be constructed from vector that does not have 9 elements"));
			}
			init_from_values( arg_vector[ 0 ], arg_vector[ 1 ], arg_vector[ 2 ],
			                  arg_vector[ 3 ], arg_vector[ 4 ], arg_vector[ 5 ],
			                  arg_vector[ 6 ], arg_vector[ 7 ], arg_vector[ 8 ] );
			return *this;
		}

		/// \brief Private initialisation from row-major order vector to be used by multiple ctors
		///        and then basic checks that this is a valid rotation
		inline rotation & rotation::init_from_values_with_rotation_checks(const double &arg_val_00, const double &arg_val_01, const double &arg_val_02,
		                                                                  const double &arg_val_10, const double &arg_val_11, const double &arg_val_12,
		                                                                  const double &arg_val_20, const double &arg_val_21, const double &arg_val_22
		                                                                  ) {
			init_from_values( arg_val_00, arg_val_01, arg_val_02,
			                  arg_val_10, arg_val_11, arg_val_12,
			                  arg_val_20, arg_val_21, arg_val_22 );
			check_is_valid_rotation();
			return *this;
		}

		/// \brief Private initialisation from row-major order vector to be used by multiple ctors
		inline rotation & rotation::init_from_values(const double &arg_val_00, const double &arg_val_01, const double &arg_val_02,
		                                             const double &arg_val_10, const double &arg_val_11, const double &arg_val_12,
		                                             const double &arg_val_20, const double &arg_val_21, const double &arg_val_22
		                                             ) {
			value_0_0 = arg_val_00;
			value_0_1 = arg_val_01;
			value_0_2 = arg_val_02;
			value_1_0 = arg_val_10;
			value_1_1 = arg_val_11;
			value_1_2 = arg_val_12;
			value_2_0 = arg_val_20;
			value_2_1 = arg_val_21;
			value_2_2 = arg_val_22;
			return *this;
		}

		/// \brief Set a specific value from a rotation object's rotation matrix
		inline rotation & rotation::set_value(const size_t &arg_row_index, ///< TODOCUMENT
		                                      const size_t &arg_col_index, ///< TODOCUMENT
		                                      const double &arg_value      ///< TODOCUMENT
		                                      ) {
			check_index( arg_row_index );
			check_index( arg_col_index );
			check_value( arg_value     );
			if      ( arg_row_index == 0 && arg_col_index == 0 ) {
				value_0_0 = arg_value;
			}
			else if ( arg_row_index == 0 && arg_col_index == 1 ) {
				value_0_1 = arg_value;
			}
			else if ( arg_row_index == 0 && arg_col_index == 2 ) {
				value_0_2 = arg_value;
			}
			else if ( arg_row_index == 1 && arg_col_index == 0 ) {
				value_1_0 = arg_value;
			}
			else if ( arg_row_index == 1 && arg_col_index == 1 ) {
				value_1_1 = arg_value;
			}
			else if ( arg_row_index == 1 && arg_col_index == 2 ) {
				value_1_2 = arg_value;
			}
			else if ( arg_row_index == 2 && arg_col_index == 0 ) {
				value_2_0 = arg_value;
			}
			else if ( arg_row_index == 2 && arg_col_index == 1 ) {
				value_2_1 = arg_value;
			}
			else if ( arg_row_index == 2 && arg_col_index == 2 ) {
				value_2_2 = arg_value;
			}
			BOOST_THROW_EXCEPTION(cath::common::out_of_range_exception("Rotation row and/or column indices unrecognised"));
			return *this;
		}

		/// \brief TODOCUMENT
		inline void rotation::check_is_valid_rotation() const {
#ifndef NDEBUG
		//	using boost::math::isfinite;
		//	for (const double_vec &row : matrix) {
		//		for (const double &element : row) {
		//			if (!boost::math::isfinite(element)) {
		//				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("All elements of the rotation must be a normal, finite floating-point numbers"));
		//			}
		//		}
		//	}
			if ( ! boost::math::isnormal( tolerance ) || tolerance >= 1.0 || tolerance <= 0.0 ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("The tolerance of a rotation must be a normal value between 0.0 and 1.0"));
			}

			// Check that the determinant is (close to) 1.0
			if ( cath::common::difference( determinant(), 1.0 ) > tolerance ) {
				BOOST_THROW_EXCEPTION(
					cath::common::runtime_error_exception(
						"Invalid rotation: determinant is " + boost::lexical_cast<std::string>( determinant() ) + ", which is not close enough to 1"
					)
				);
			}

			// Check that the multiplying by the transpose gives (close to) the identity matrix
			//
			// NOTE: Cannot just call transpose_of() here because that uses a constructor which starts an infinite recursion
			for (const size_t &check_row_ctr : common::indices( coord::NUM_DIMS ) ) {
				for (const size_t &check_col_ctr : common::indices( coord::NUM_DIMS ) ) {
					// Can't use IDENTITY_ROTATION here because these call may have been triggered by IDENTITY_ROTATION's ctor
					const double correct_value = (check_row_ctr == check_col_ctr) ? 1.0 : 0.0;

					const double transpose_row_mult_value = dot_product(
						get_row( *this, check_row_ctr ),
						get_row( *this, check_col_ctr )
					);
					const double transpose_row_mult_error = cath::common::difference( transpose_row_mult_value, correct_value );
					const double transpose_col_mult_value = dot_product(
						get_column( *this, check_row_ctr ),
						get_column( *this, check_col_ctr )
					);
					const double transpose_col_mult_error = cath::common::difference( transpose_col_mult_value, correct_value );


					if (transpose_row_mult_error > tolerance) {
						BOOST_THROW_EXCEPTION(cath::common::runtime_error_exception(
							"Invalid rotation: when its transpose is multiplied by it, the row result is not close enough to the identity matrix (error is "
							+ boost::lexical_cast<std::string>(transpose_row_mult_error)
							+ " and tolerance is "
							+ boost::lexical_cast<std::string>(tolerance)
							+ ")"
						));
					}
					if (transpose_col_mult_error > tolerance) {
						BOOST_THROW_EXCEPTION(cath::common::runtime_error_exception(
							"Invalid rotation: when its transpose is multiplied by it, the column result is not close enough to the identity matrix (error is "
							+ boost::lexical_cast<std::string>(transpose_col_mult_error)
							+ " and tolerance is "
							+ boost::lexical_cast<std::string>(tolerance)
							+ ")"
						));
					}
				}
			}
#endif
		}

		/// \brief TODOCUMENT
		inline double rotation::determinant() const {
			return 0.0
				+ get_value<0, 0>() * ( get_value<1, 1>() * get_value<2, 2>() - get_value<1, 2>() * get_value<2, 1>() )
				+ get_value<0, 1>() * ( get_value<1, 2>() * get_value<2, 0>() - get_value<1, 0>() * get_value<2, 2>() )
				+ get_value<0, 2>() * ( get_value<1, 0>() * get_value<2, 1>() - get_value<1, 1>() * get_value<2, 0>() );
		}

		/// \brief TODOCUMENT
		template <size_t index>
		inline void rotation::check_index() const {
			static_assert( index >= 0 && index <= 2, "Rotation index template parameter is not between 0 and 2 inclusive" );
		}

		/// \brief TODOCUMENT
		template <size_t row_index, size_t col_index>
		inline void rotation::set_value(const double &arg_value ///< TODOCUMENT
		                                ) {
			check_index<row_index>();
			check_index<col_index>();
			if      ( row_index == 0 && col_index == 0 ) {
				value_0_0 = arg_value;
			}
			else if ( row_index == 0 && col_index == 1 ) {
				value_0_1 = arg_value;
			}
			else if ( row_index == 0 && col_index == 2 ) {
				value_0_2 = arg_value;
			}
			else if ( row_index == 1 && col_index == 0 ) {
				value_1_0 = arg_value;
			}
			else if ( row_index == 1 && col_index == 1 ) {
				value_1_1 = arg_value;
			}
			else if ( row_index == 1 && col_index == 2 ) {
				value_1_2 = arg_value;
			}
			else if ( row_index == 2 && col_index == 0 ) {
				value_2_0 = arg_value;
			}
			else if ( row_index == 2 && col_index == 1 ) {
				value_2_1 = arg_value;
			}
			else if ( row_index == 2 && col_index == 2 ) {
				value_2_2 = arg_value;
			}
			BOOST_THROW_EXCEPTION(cath::common::out_of_range_exception("Rotation row and/or column indices unrecognised"));
		}

		/// \brief Constructor for rotation using row-major order
		inline rotation::rotation(const double &val_00, const double &val_01, const double &val_02,
		                          const double &val_10, const double &val_11, const double &val_12,
		                          const double &val_20, const double &val_21, const double &val_22,
		                          const double &arg_tolerance
		                          ) : tolerance( arg_tolerance ) {
			init_from_values_with_rotation_checks(val_00, val_01, val_02,
			                                      val_10, val_11, val_12,
			                                      val_20, val_21, val_22 );
		}

		/// \brief Constructor from row-major order vector
		inline rotation::rotation(const doub_vec &arg_vector, ///< TODOCUMENT
		                          const double   &arg_tolerance
		                          ) : tolerance(arg_tolerance) {
			init_from_vector_with_rotation_checks( arg_vector );
		}

		/// \brief Get a specific value from a rotation object's rotation matrix
		inline const double & rotation::get_value(const size_t &arg_row_index, ///< TODOCUMENT
		                                          const size_t &arg_col_index  ///< TODOCUMENT
		                                          ) const {
			check_index( arg_row_index );
			check_index( arg_col_index );
			if ( arg_row_index == 0 && arg_col_index == 0 ) {
				return value_0_0;
			}
			if ( arg_row_index == 0 && arg_col_index == 1 ) {
				return value_0_1;
			}
			if ( arg_row_index == 0 && arg_col_index == 2 ) {
				return value_0_2;
			}
			if ( arg_row_index == 1 && arg_col_index == 0 ) {
				return value_1_0;
			}
			if ( arg_row_index == 1 && arg_col_index == 1 ) {
				return value_1_1;
			}
			if ( arg_row_index == 1 && arg_col_index == 2 ) {
				return value_1_2;
			}
			if ( arg_row_index == 2 && arg_col_index == 0 ) {
				return value_2_0;
			}
			if ( arg_row_index == 2 && arg_col_index == 1 ) {
				return value_2_1;
			}
			if ( arg_row_index == 2 && arg_col_index == 2 ) {
				return value_2_2;
			}
			BOOST_THROW_EXCEPTION(cath::common::out_of_range_exception("Rotation row and/or column indices unrecognised"));
			return value_0_0; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}

		/// \brief TODOCUMENT
		template <size_t row_index, size_t col_index>
		inline const double & rotation::get_value() const {
			check_index<row_index>();
			check_index<col_index>();
			if ( row_index == 0 && col_index == 0 ) {
				return value_0_0;
			}
			if ( row_index == 0 && col_index == 1 ) {
				return value_0_1;
			}
			if ( row_index == 0 && col_index == 2 ) {
				return value_0_2;
			}
			if ( row_index == 1 && col_index == 0 ) {
				return value_1_0;
			}
			if ( row_index == 1 && col_index == 1 ) {
				return value_1_1;
			}
			if ( row_index == 1 && col_index == 2 ) {
				return value_1_2;
			}
			if ( row_index == 2 && col_index == 0 ) {
				return value_2_0;
			}
			if ( row_index == 2 && col_index == 1 ) {
				return value_2_1;
			}
			if ( row_index == 2 && col_index == 2 ) {
				return value_2_2;
			}
		}


		/// \brief Multiply this rotation on the right by arg_rotation (ie (*this) = (*this) * arg_rotation )
		///
		/// boost::multipliable<rotation> also uses this to add the non-member function:
		///  rotation operator*(const rotation &, const rotation &);
		///
		/// Remember that rotations are associative (so it is generally true that A*(B*C) = (A*B)*C)
		/// but not commutative (so it is not generally true that A*B = B*A)
		inline void rotation::operator*=(const rotation &arg_rotation ///< TODOCUMENT
		                                 ) {
			const double val_00 = dot_product( get_row<0>( *this ), get_column<0>( arg_rotation ) );
			const double val_01 = dot_product( get_row<0>( *this ), get_column<1>( arg_rotation ) );
			const double val_02 = dot_product( get_row<0>( *this ), get_column<2>( arg_rotation ) );

			const double val_10 = dot_product( get_row<1>( *this ), get_column<0>( arg_rotation ) );
			const double val_11 = dot_product( get_row<1>( *this ), get_column<1>( arg_rotation ) );
			const double val_12 = dot_product( get_row<1>( *this ), get_column<2>( arg_rotation ) );

			const double val_20 = dot_product( get_row<2>( *this ), get_column<0>( arg_rotation ) );
			const double val_21 = dot_product( get_row<2>( *this ), get_column<1>( arg_rotation ) );
			const double val_22 = dot_product( get_row<2>( *this ), get_column<2>( arg_rotation ) );

			init_from_values( val_00, val_01, val_02,
			                  val_10, val_11, val_12,
			                  val_20, val_21, val_22 );
		}

		/// \brief Get a specific row from a rotation object's rotation matrix
		///
		/// \relates rotation
		template<size_t row_index>
		inline coord get_row(const rotation &arg_rotation ///< TODOCUMENT
		                     ) {
			return coord(
				arg_rotation.get_value<row_index, 0>(),
				arg_rotation.get_value<row_index, 1>(),
				arg_rotation.get_value<row_index, 2>()
			);
		}

		/// \brief Get a specific column from a rotation object's rotation matrix
		///
		/// \relates rotation
		template<size_t col_index>
		inline coord get_column(const rotation &arg_rotation ///< TODOCUMENT
		                        ) {
			return coord(
				arg_rotation.get_value<0, col_index>(),
				arg_rotation.get_value<1, col_index>(),
				arg_rotation.get_value<2, col_index>()
			);
		}

		/// \brief Get a specific row from a rotation object's rotation matrix
		///
		/// \relates rotation
		inline coord get_row(const rotation &arg_rotation, ///< TODOCUMENT
		                     const size_t   &arg_row_index ///< TODOCUMENT
		                     ) {
			return coord(
				arg_rotation.get_value( arg_row_index, 0 ),
				arg_rotation.get_value( arg_row_index, 1 ),
				arg_rotation.get_value( arg_row_index, 2 )
			);
		}

		/// \brief Get a specific column from a rotation object's rotation matrix
		///
		/// \relates rotation
		inline coord get_column(const rotation &arg_rotation,    ///< TODOCUMENT
		                        const size_t   &arg_column_index ///< TODOCUMENT
		                        ) {
			return coord(
				arg_rotation.get_value( 0, arg_column_index ),
				arg_rotation.get_value( 1, arg_column_index ),
				arg_rotation.get_value( 2, arg_column_index )
			);
		}

		/// \brief Return the result of transposing this rotation.
		///
		/// \relates rotation
		///
		/// Since the input should be constrained to be a valid rotation, this will also act as the inverse rotation
		inline rotation transpose_copy(const rotation &arg_rotation ///< TODOCUMENT
		                               ) {
			return rotation(
				arg_rotation.get_value<0, 0>(), arg_rotation.get_value<1, 0>(), arg_rotation.get_value<2, 0>(),
				arg_rotation.get_value<0, 1>(), arg_rotation.get_value<1, 1>(), arg_rotation.get_value<2, 1>(),
				arg_rotation.get_value<0, 2>(), arg_rotation.get_value<1, 2>(), arg_rotation.get_value<2, 2>()
			);
		}

		/// \brief Calculate the rotation that gives the second rotation if applied after the first rotation
		///
		/// \relates rotation
		inline rotation rotation_between_rotations(const rotation &arg_rotation_1, ///< The first rotation
		                                           const rotation &arg_rotation_2  ///< The second rotation
		                                           ) {
			return arg_rotation_2 * transpose_copy(arg_rotation_1);
		}

		/// \brief TODOCUMENT
		///
		/// \relates rotation
		inline double trace(const rotation &arg_rotation ///< TODOCUMENT
		                    ) {
			return arg_rotation.get_value<0, 0>()
			     + arg_rotation.get_value<1, 1>()
				 + arg_rotation.get_value<2, 2>();
		}

		/// \brief Calculate the angle (in radians) of the rotation
		///
		/// \relates rotation
		inline doub_angle angle_of_rotation(const rotation &arg_rotation ///< TODOCUMENT
		                                    ) {
			const double rotation_trace = trace( arg_rotation );
//			std::cerr << "rotation_trace " << rotation_trace << std::endl;

			return make_angle_from_radians<double>( acos( boost::algorithm::clamp( ( rotation_trace - 1.0 ) / 2.0, -1.0, 1.0 ) ) );
		}

		/// \brief Calculate the angle (in radians) between the two rotations
		///
		/// \relates rotation
		inline doub_angle angle_between_rotations(const rotation &arg_rotation_1, ///< The first rotation
		                                          const rotation &arg_rotation_2  ///< The second rotation
		                                          ) {
//			std::cerr << "arg_rotation_1 : " << arg_rotation_1 << std::endl;
//			std::cerr << "arg_rotation_2 : " << arg_rotation_2 << std::endl;
			return angle_of_rotation( rotation_between_rotations( arg_rotation_1, arg_rotation_2 ) );
		}


		/// \brief Rotate a coord_list object based on the specification in a rotation object
		///
		/// \relates rotation
		///
		/// \todo I think this is a good candidate for making into an operator*() non-member function
		///       Will need to be careful with symmetry if using boost::operators to provide it via a member function.
		inline void rotate(const rotation &arg_rotation,  ///< TODOCUMENT
		                   coord_list     &arg_coord_list ///< TODOCUMENT
		                   ) {
			for (coord &loop_coord : arg_coord_list) {
				::cath::geom::rotate( arg_rotation, loop_coord );
			}
		}

		/// \brief Rotate a coord object based on the specification in a rotation object
		///
		/// \relates rotation
		///
		/// \todo I think this is a good candidate for making into an operator*() non-member function
		///       Will need to be careful with symmetry if using boost::operators to provide it via a member function.
		inline void rotate(const rotation &arg_rotation, ///< TODOCUMENT
		                   coord          &arg_coord     ///< TODOCUMENT
		                   ) {
			arg_coord = coord(
				dot_product( get_row<0>( arg_rotation ), arg_coord ),
				dot_product( get_row<1>( arg_rotation ), arg_coord ),
				dot_product( get_row<2>( arg_rotation ), arg_coord )
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates rotation
		inline coord_list rotate_copy(const rotation   &arg_rotation,  ///< TODOCUMENT
		                              const coord_list &arg_coord_list ///< TODOCUMENT
		                              ) {
			coord_list new_coord_list( arg_coord_list );
			::cath::geom::rotate( arg_rotation, new_coord_list );
			return new_coord_list;
		}

		/// \brief TODOCUMENT
		///
		/// \relates rotation
		inline coord rotate_copy(const rotation &arg_rotation, ///< TODOCUMENT
		                         const coord    &arg_coord     ///< TODOCUMENT
		                         ) {
			coord new_coord( arg_coord );
			::cath::geom::rotate( arg_rotation, new_coord );
			return new_coord;
		}

	} // namespace geom

	namespace common {

		/// \brief Specialisation of cath::common::read_from_ptree for rotation
		template <>
		inline geom::rotation read_from_ptree<geom::rotation>(const boost::property_tree::ptree &arg_ptree ///< The ptree from which to read the rotation
		                                                      ) {
			return geom::rotation_from_ptree( arg_ptree );
		}

	} // namespace common
} // namespace cath

#endif
