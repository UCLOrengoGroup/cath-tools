/// \file
/// \brief The quat_rot class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_QUAT_ROT_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_QUAT_ROT_HPP

#include <boost/math/quaternion.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/mpl/vector.hpp>

#include "cath/structure/geometry/rotation.hpp"

#include <algorithm>
#include <cmath>
#include <iostream> // ***** TEMPORARY *****
#include <random>

namespace cath::geom {

	/// \brief TODOCUMENT
	///
	/// Quaternion-based rotation
	///
	/// These must always have length 1
	///
	/// \todo For when NDEBUG isn't defined, add more checking that the these quaternions
	///       all have length 1
	template <typename T>
	class quat_rot_impl final : public boost::math::quaternion<T>,
	                            private boost::multipliable<quat_rot_impl<T>, T,
	                                    boost::addable<quat_rot_impl<T>>> {
	private:
		using super = boost::math::quaternion<T>;

	public:
		explicit quat_rot_impl(const T &,
		                       const T &,
		                       const T &,
		                       const T &);
		explicit quat_rot_impl(const boost::math::quaternion<T> &);
	};

	/// \brief Ctor for quat_rot_impl
	template <typename T>
	inline quat_rot_impl<T>::quat_rot_impl(const T &prm_0, ///< TODOCUMENT
	                                       const T &prm_1, ///< TODOCUMENT
	                                       const T &prm_2, ///< TODOCUMENT
	                                       const T &prm_3  ///< TODOCUMENT
	                                       ) : super( prm_0, prm_1, prm_2, prm_3 ) {
	}

	/// \brief Ctor for quat_rot_impl
	template <typename T>
	inline quat_rot_impl<T>::quat_rot_impl(const boost::math::quaternion<T> &prm_quaternion ///< The quaternion from which this quat_rot should be constructed
	                                       ) : super( prm_quaternion / abs( prm_quaternion ) ) {
	}

	///\brief TODOCUMENT
	template <typename T>
	quat_rot_impl<T> operator-(const quat_rot_impl<T> &prm_quaternion ///< TODOCUMENT
	                           ) {
		return quat_rot_impl<T>{ quat_rot_impl<T>( 0.0, 0.0, 0.0, 0.0 ) - prm_quaternion };
	}

	/// \brief Make the identity quat_rot
	template <typename T>
	inline quat_rot_impl<T> make_identity_quat_rot() {
		return quat_rot_impl<T>( boost::math::quaternion<T>( 1.0, 0.0, 0.0, 0.0 ) );
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline void normalise(quat_rot_impl<T> &prm_quaternion ///< TODOCUMENT
	                      ) {
		prm_quaternion /= abs( prm_quaternion );
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline quat_rot_impl<T> normalise_copy(quat_rot_impl<T> prm_quaternion ///< TODOCUMENT
	                                       ) {
		normalise( prm_quaternion );
		return prm_quaternion;
	}

	// /// \brief TODOCUMENT
	// template <typename T>
	// inline quat_rot_impl<T> conj(quat_rot_impl<T> prm_quaternion ///< TODOCUMENT
	//                              ) {
	// 	return conj( prm_quaternion );
	// }

	/// \brief TODOCUMENT
	template <typename T, typename U>
	inline void rotate(const quat_rot_impl<T>                                              &prm_quaternion, ///< TODOCUMENT
	                   boost::geometry::model::point<U, 3, boost::geometry::cs::cartesian> &prm_point       ///< TODOCUMENT
	                   ) {
		const auto point = quat_rot_impl<T>(
			0.0,
			boost::geometry::get<0>( prm_point ),
			boost::geometry::get<1>( prm_point ),
			boost::geometry::get<2>( prm_point )
		);
		const auto rot_quat = prm_quaternion * point * conj( prm_quaternion );
		boost::geometry::set<0>( prm_point, rot_quat.R_component_2() );
		boost::geometry::set<1>( prm_point, rot_quat.R_component_3() );
		boost::geometry::set<2>( prm_point, rot_quat.R_component_4() );
	}

	/// \brief TODOCUMENT
	template <typename T, typename U>
	inline boost::geometry::model::point<U, 3, boost::geometry::cs::cartesian> rotate_copy(const quat_rot_impl<T>                                              &prm_quaternion, ///< TODOCUMENT
	                                                                                       boost::geometry::model::point<U, 3, boost::geometry::cs::cartesian>  prm_point       ///< TODOCUMENT
	                                                                                       ) {
		::cath::geom::rotate( prm_quaternion, prm_point );
		return prm_point;
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline T inner_product(const quat_rot_impl<T> &prm_quat_rot_a, ///< TODOCUMENT
	                       const quat_rot_impl<T> &prm_quat_rot_b  ///< TODOCUMENT
	                       ) {
		return (
			  prm_quat_rot_a.R_component_1() * prm_quat_rot_b.R_component_1()
			+ prm_quat_rot_a.R_component_2() * prm_quat_rot_b.R_component_2()
			+ prm_quat_rot_a.R_component_3() * prm_quat_rot_b.R_component_3()
			+ prm_quat_rot_a.R_component_4() * prm_quat_rot_b.R_component_4()
		);
	}

	/// \brief Make a quat_rot that represents the same rotation as the specified rotation
	template <typename T>
	inline quat_rot_impl<T> make_quat_rot_from_rotation(const rotation &prm_rotation ///< The matrix rotation from which the quat_rot should be built
	                                                    ) {
		const T  zero = static_cast<T>( 0.0 );
		const T  one  = static_cast<T>( 1.0 );
		const T  half = static_cast<T>( 0.5 );

		const T &q_xx = static_cast<T>( prm_rotation.get_value<0, 0>() );
		const T &q_xy = static_cast<T>( prm_rotation.get_value<0, 1>() );
		const T &q_xz = static_cast<T>( prm_rotation.get_value<0, 2>() );

		const T &q_yx = static_cast<T>( prm_rotation.get_value<1, 0>() );
		const T &q_yy = static_cast<T>( prm_rotation.get_value<1, 1>() );
		const T &q_yz = static_cast<T>( prm_rotation.get_value<1, 2>() );

		const T &q_zx = static_cast<T>( prm_rotation.get_value<2, 0>() );
		const T &q_zy = static_cast<T>( prm_rotation.get_value<2, 1>() );
		const T &q_zz = static_cast<T>( prm_rotation.get_value<2, 2>() );

		// const T  t    = q_xx + q_yy + q_zz;
		// const T  r    = std::sqrt( one + t );
		const T  w    =                half * std::sqrt( std::max( zero, one + q_xx + q_yy + q_zz ) );
		const T  x    = std::copysign( half * std::sqrt( std::max( zero, one + q_xx - q_yy - q_zz ) ), q_zy - q_yz );
		const T  y    = std::copysign( half * std::sqrt( std::max( zero, one - q_xx + q_yy - q_zz ) ), q_xz - q_zx );
		const T  z    = std::copysign( half * std::sqrt( std::max( zero, one - q_xx - q_yy + q_zz ) ), q_yx - q_xy );

		// std::cerr << "q_xx is : " << q_xx << std::endl;
		// std::cerr << "q_xy is : " << q_xy << std::endl;
		// std::cerr << "q_xz is : " << q_xz << std::endl;
		// std::cerr << "t is    : " << t    << std::endl;
		// std::cerr << "r is    : " << r    << std::endl;
		// std::cerr << "w is    : " << w    << std::endl;
		// std::cerr << "x is    : " << x    << std::endl;
		// std::cerr << "y is    : " << y    << std::endl;
		// std::cerr << "z is    : " << z    << std::endl;

		return quat_rot_impl<T>( boost::math::quaternion<T>( w, x, y, z ) );
	}

	/// \brief Make a rotation that represents the same rotation as the specified quat_rot
	template <typename T>
	inline rotation make_rotation_from_quat_rot(const quat_rot_impl<T> &prm_quat_rot ///< The quat_rot from which the rotation should be built
	                                            ) {
		const T one = static_cast<T>( 1.0 );
		const T two = static_cast<T>( 2.0 );
		const T q0  = prm_quat_rot.R_component_1();
		const T q1  = prm_quat_rot.R_component_2();
		const T q2  = prm_quat_rot.R_component_3();
		const T q3  = prm_quat_rot.R_component_4();
		return rotation(
			boost::numeric_cast<double>( one - two * q2 * q2 - two * q3 * q3 ),
			boost::numeric_cast<double>(       two * q1 * q2 - two * q0 * q3 ),
			boost::numeric_cast<double>(       two * q1 * q3 + two * q0 * q2 ),

			boost::numeric_cast<double>(       two * q2 * q1 + two * q0 * q3 ),
			boost::numeric_cast<double>( one - two * q3 * q3 - two * q1 * q1 ),
			boost::numeric_cast<double>(       two * q2 * q3 - two * q0 * q1 ),

			boost::numeric_cast<double>(       two * q3 * q1 - two * q0 * q2 ),
			boost::numeric_cast<double>(       two * q3 * q2 + two * q0 * q1 ),
			boost::numeric_cast<double>( one - two * q1 * q1 - two * q2 * q2 )
		);
	}

	/// \brief Calculate the angle required to perform the rotation represented by the specified quat_rot
	template <typename T>
	inline angle<T> angle_of_quat_rot(const quat_rot_impl<T> &prm_quat_rot ///< The quat_rot
	                                  ) {
		return make_angle_from_radians<T>(
			boost::numeric_cast<double>(
				2.0 * std::acos( std::min(
					static_cast<T>( 1.0 ),
					std::fabs( prm_quat_rot.R_component_1() )
				) )
			)
		);
	}

	/// \brief Calculate the angle required to rotate between the two orientations represented by the specified quat_rots
	template <typename T>
	inline angle<T> angle_between_quat_rots(const quat_rot_impl<T> &prm_quat_rot_a, ///< The first quat_rot
	                                        const quat_rot_impl<T> &prm_quat_rot_b  ///< The second quat_rot
	                                        ) {
		return angle_of_quat_rot<T>( quat_rot_impl<T>( prm_quat_rot_b * conj( prm_quat_rot_a ) ) );
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline quat_rot_impl<T> rotation_between_rotations(const quat_rot_impl<T> &prm_quat_rot_a, ///< The first quat_rot
	                                                   const quat_rot_impl<T> &prm_quat_rot_b  ///< The second quat_rot
	                                                   ) {
		// return quat_rot_impl<T>{ conj( prm_quat_rot_a ) *       prm_quat_rot_b   };
		// return quat_rot_impl<T>{       prm_quat_rot_a   * conj( prm_quat_rot_b ) };
		// return quat_rot_impl<T>{ conj( prm_quat_rot_b ) *       prm_quat_rot_a   };
		return quat_rot_impl<T>{       prm_quat_rot_b   * conj( prm_quat_rot_a ) };
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline quat_rot_impl<T> interpolate_smooth_line(const quat_rot_impl<T> &prm_quat_rot_a, ///< The first quat_rot
	                                                const quat_rot_impl<T> &prm_quat_rot_b, ///< The second quat_rot
	                                                const T                &prm_frac_from_a ///< The (linear) fraction from a to b of the required quaternion
	                                                ) {
		if ( prm_frac_from_a < 0.0 || prm_frac_from_a > 1.0 ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot interpolate quaternions with value outside range [0,1]"));
		}
		const auto dot            = inner_product( prm_quat_rot_a, prm_quat_rot_b );
		const auto frac_remaining = static_cast<T>( 1.0 ) - prm_frac_from_a;
		return ( dot < 0.0 )
			? normalise_copy( ( frac_remaining * prm_quat_rot_a ) + ( prm_frac_from_a * ( - prm_quat_rot_b ) ) )
			: normalise_copy( ( frac_remaining * prm_quat_rot_a ) + ( prm_frac_from_a *     prm_quat_rot_b   ) );
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline quat_rot_impl<T> interpolate_smooth_angle(const quat_rot_impl<T> &prm_quat_rot_a, ///< The first quat_rot
	                                                 const quat_rot_impl<T> &prm_quat_rot_b, ///< The second quat_rot
	                                                 const T                &prm_frac_from_a ///< The (angular) fraction from a to b of the required quaternion
	                                                 ) {
		const auto frac_remaining         = static_cast<T>( 1.0 ) - prm_frac_from_a;
		const auto half_angle_between     = 0.5 * angle_between_quat_rots<T>( prm_quat_rot_a, prm_quat_rot_b );
		const auto this_part              = std::sin( angle_in_radians( prm_frac_from_a * half_angle_between ) );
		const auto that_part              = std::sin( angle_in_radians( frac_remaining  * half_angle_between ) );
		return interpolate_smooth_line(
			prm_quat_rot_a,
			prm_quat_rot_b,
			this_part / ( this_part + that_part )
		);
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline quat_rot_impl<T> mid_point(const quat_rot_impl<T> &prm_quat_rot_a, ///< The first  quat_rot
	                                  const quat_rot_impl<T> &prm_quat_rot_b  ///< The second quat_rot
	                                  ) {
		return interpolate_smooth_line( prm_quat_rot_a, prm_quat_rot_b, static_cast<T>( 0.5 ) );
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline quat_rot_impl<T> from_first_toward_second_at_angle(const quat_rot_impl<T> &prm_quat_rot_a,     ///< The first  quat_rot
	                                                          const quat_rot_impl<T> &prm_quat_rot_b,     ///< The second quat_rot
	                                                          const angle<T>         &prm_requested_angle ///< TODOCUMENT
	                                                          ) {
		const auto angle_between = angle_between_quat_rots<T>( prm_quat_rot_a, prm_quat_rot_b );
		if ( prm_requested_angle > angle_between ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to generate quat_rot at angle x from a to b if x > the angle between a and b"));
		}
		if ( angle_between.get_angle() == 0.0 ) {
			return prm_quat_rot_a;
		}
		const T ratio = prm_requested_angle / angle_between;
		return interpolate_smooth_angle( prm_quat_rot_a, prm_quat_rot_b, ratio );
	}

	/// \brief Calculate the distance between the two orientations represented by the specified quat_rots
	///        ( 1 - a.b )
	///
	/// This should be a bit faster than the full calculation of the angle between quat_rots
	template <typename T>
	inline T distance_1_between_quat_rots(const quat_rot_impl<T> &prm_quat_rot_a, ///< The first quat_rot
	                                      const quat_rot_impl<T> &prm_quat_rot_b  ///< The second quat_rot
	                                      ) {
		return static_cast<T>( 1.0 ) - std::fabs( inner_product( prm_quat_rot_a, prm_quat_rot_b ) );
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline T distance_1_of_quat_rot(const quat_rot_impl<T> &prm_quat_rot ///< The quat_rot
	                                ) {
		return distance_1_between_quat_rots<T>(
			make_identity_quat_rot<T>(),
			prm_quat_rot
		);
	}

	/// \brief TODOCUMENT
	template <typename T, typename U>
	inline T distance_1_of_angle(const angle<U> &prm_angle ///< TODOCUMENT
	                             ) {
		return distance_1_of_quat_rot<T>(
			make_quat_rot_from_rotation<T>(
				rotation_of_angle( prm_angle )
			)
		);
	}

	/// \brief Generate a random quat_rot (which must of course have length 1)
	template <typename T>
	quat_rot_impl<T> make_random_quat_rot(std::default_random_engine &prm_rng ///< TODOCUMENT
	                                      ) {
		// Create a uniform distribution for generating the four parts of the quaternion
		std::uniform_real_distribution<T> the_dist( -1.0, 1.0 );

		// Construct the parts and make a
		// (note: the random parts are generated first for the sake of reproducibility
		//  because generating them as arguments to the quat_rot ctor would mean that the
		//  compiler was free to arbitrarily reorder their generations)
		const T q0 = the_dist( prm_rng );
		const T q1 = the_dist( prm_rng );
		const T q2 = the_dist( prm_rng );
		const T q3 = the_dist( prm_rng );
		const auto raw_quat_rot = quat_rot_impl<T>{ q0, q1, q2, q3 };

		// In the very unlikely circumstance that all parts are 0.0, recurse
		// to try again
		if ( q0 == 0.0 && q1 == 0.0 && q2 == 0.0 && q3 == 0.0 ) {
			return make_random_quat_rot<T>( prm_rng );
		}

		// Return a normalised copy of the raw quaternion
		return normalise_copy( raw_quat_rot );
	}

	/// \brief TODOCUMENT
	template <typename T>
	std::ostream & operator<<(std::ostream           &prm_os,      ///< TODOCUMENT
	                          const quat_rot_impl<T> &prm_quat_rot ///< TODOCUMENT
	                          ) {
		std::ostringstream temp_ss;
		temp_ss << "quat_rot[";
		temp_ss << std::right << std::setw( 7 ) << prm_quat_rot.R_component_1();
		temp_ss << ",";
		temp_ss << std::right << std::setw( 7 ) << prm_quat_rot.R_component_2();
		temp_ss << ",";
		temp_ss << std::right << std::setw( 7 ) << prm_quat_rot.R_component_3();
		temp_ss << ",";
		temp_ss << std::right << std::setw( 7 ) << prm_quat_rot.R_component_4();
		temp_ss << "]";
		prm_os << temp_ss.str();
		return prm_os;
	}

	// /// \brief TODOCUMENT
	// template <typename T>
	// inline void rotate(const quat_rot_impl<T> &prm_quat_rot, ///< TODOCUMENT
	//                    coord                  &prm_coord     ///< TODOCUMENT
	//                    ) {
	// 	prm_quat_rot;
	// 	prm_coord;
	// }

	// /// \brief TODOCUMENT
	// template <typename T>
	// inline void rotate(const quat_rot_impl<T> &prm_quat_rot,  ///< TODOCUMENT
	//                    coord_list             &prm_coord_list ///< TODOCUMENT
	//                    ) {
	// 	for (const coord &the_coord : prm_coord_list) {
	// 		rotate( prm_quat_rot, the_coord );
	// 	}
	// }

	// /// \brief TODOCUMENT
	// template <typename T>
	// inline coord rotate_copy(const quat_rot_impl<T> &prm_quat_rot, ///< TODOCUMENT
	//                          coord                   prm_coord    ///< TODOCUMENT
	//                          ) {
	// 	rotate( prm_quat_rot, prm_coord );
	// 	return prm_coord;
	// }

	// /// \brief TODOCUMENT
	// template <typename T>
	// inline coord_list rotate_copy(const quat_rot_impl<T> &prm_quat_rot,  ///< TODOCUMENT
	//                               coord_list              prm_coord_list ///< TODOCUMENT
	//                               ) {
	// 	rotate( prm_quat_rot, prm_coord_list );
	// 	return prm_coord_list;
	// }

	using all_quat_rot_types = boost::mpl::vector<float>;
	// using all_quat_rot_types = boost::mpl::vector<       double, long double>;
	// using quat_rot = quat_rot_impl<double>;
	using float_quat_rot = quat_rot_impl<float>;
} // namespace cath::geom

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_QUAT_ROT_HPP
