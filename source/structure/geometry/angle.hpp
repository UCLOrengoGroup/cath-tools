/// \file
/// \brief The angle class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_ANGLE_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_ANGLE_H

#include <boost/math/constants/constants.hpp>
#include <boost/operators.hpp>
//#include <boost/units/quantity.hpp>
//#include <boost/units/systems/si/plane_angle.hpp>

#include "common/difference.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "structure/geometry/angle_endpoint_loc.hpp"

#include <algorithm>
#include <cmath>
// #include <iostream>

namespace cath {
	namespace geom {
		template <typename T>
		class angle;

		template<typename T>
		inline angle<T> make_angle_from_degrees(const double &);

		/// \brief TODOCUMENT
		template <typename T>
		class angle final : private boost::additive       <angle<T>,
		                            boost::totally_ordered<angle<T>,
		                            boost::multiplicative <angle<T>, T> > > {
		private:
			using angle_type = T;
//			using angle_type = boost::units::quantity<boost::units::si::plane_angle>;

			/// \brief TODOCUMENT
			angle_type angle_value;

			angle_type & get_angle();

		public:
			explicit angle(const angle_type &);

			const angle_type & get_angle() const;

			angle & operator-=(const angle &);
			angle & operator+=(const angle &);
			angle & operator*=(const T &);
			angle & operator/=(const T &);

			angle & quick_shift();

			angle & shift(const angle & = make_angle_from_degrees<angle_type>( 0.0 ),
			              const angle_endpoint_loc & = angle_endpoint_loc::USE_LOWER);
		};

		template <typename T>
		angle<T> operator-(const angle<T> &);

		template <typename T>
		T operator/(const angle<T> &,
		            const angle<T> &);

		template <typename T>
		bool operator<(const angle<T> &,
		               const angle<T> &);

		template <typename T>
		bool operator==(const angle<T> &,
		                const angle<T> &);

		template <typename T>
		std::ostream & operator<<(std::ostream &,
		                          const angle<T> &);

		template <typename T>
		inline angle<T> unshifted_wrapped_difference(const angle<T> &,
		                                             const angle<T> &);

		template <typename T>
		inline angle<T> wrapped_difference(angle<T>,
		                                   angle<T>);

		template <typename T>
		angle<T> shift_copy(angle<T>,
		                    const angle<T> & = make_angle_from_degrees<T>( 0.0 ),
		                    const angle_endpoint_loc & = angle_endpoint_loc::USE_LOWER);

		template <typename T>
		inline angle<T> make_angle_from_degrees(const double &);
		template <typename T>
		inline angle<T> make_angle_from_radians(const double &);
		template <typename T>
		inline angle<T> make_angle_from_revolutions(const double &);

		template <typename T>
		T angle_in_degrees(const angle<T> &);

		template <typename T>
		T angle_in_radians(const angle<T> &);

		template <typename T>
		T angle_in_revolutions(const angle<T> &);

		template <typename T>
		angle<T> zero_angle();
		template <typename T>
		angle<T> one_revolution();


		/// \brief Ctor for angle
		template <typename T>
		inline angle<T>::angle(const angle_type &arg_angle_value ///< TODOCUMENT
		                       ) : angle_value( arg_angle_value ) {
			if ( ! ::boost::math::isfinite( arg_angle_value ) ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Angle value must be finite"));
			}
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline typename angle<T>::angle_type & angle<T>::get_angle() {
			return angle_value;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline const typename angle<T>::angle_type & angle<T>::get_angle() const {
			return angle_value;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline angle<T> & angle<T>::operator-=(const angle<T> &arg_angle ///< TODOCUMENT
		                                       ) {
			get_angle() -= arg_angle.get_angle();
			return *this;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline angle<T> & angle<T>::operator+=(const angle<T> &arg_angle ///< TODOCUMENT
		                                       ) {
			get_angle() += arg_angle.get_angle();
			return *this;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline angle<T> & angle<T>::operator*=(const T &arg_angle ///< TODOCUMENT
		                                       ) {
			get_angle() *= arg_angle;
			return *this;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline angle<T> & angle<T>::operator/=(const T &arg_angle ///< TODOCUMENT
		                                       ) {
			get_angle() /= arg_angle;
			return *this;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline angle<T> & angle<T>::quick_shift() {
			constexpr T two_pi = static_cast<T>( 2.0 ) * boost::math::constants::pi<T>();
			while ( get_angle() < 0.0 ) {
				get_angle() += two_pi;
			}
			while ( get_angle() >= two_pi ) {
				get_angle() -= two_pi;
			}
			return *this;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline angle<T> & angle<T>::shift(const angle<T>           &arg_lower_angle_of_range, ///< TODOCUMENT
		                                  const angle_endpoint_loc &arg_endpoint_loc          ///< TODOCUMENT
		                                  ) {
			// Grab some basics
			const angle  endpoint_tolerance         = make_angle_from_degrees<T>( static_cast<T>( 0.0000000001 ) );
			const angle &lower                      = arg_lower_angle_of_range;
			const angle  upper                      = arg_lower_angle_of_range + one_revolution<T>();
			const angle  lower_less_one_rev         = lower - one_revolution<T>();
			const angle  upper_plus_one_rev         = upper + one_revolution<T>();
			const T      lower_less_one_rev_in_rads = angle_in_radians( lower_less_one_rev );
			const T      upper_plus_one_rev_in_rads = angle_in_radians( upper_plus_one_rev );

			// Get the angle into the rough area (to do most of the work if it's very far away)
			// (explanation:
			//   if more than one revolution above required range, put it into the range that's [0, 2pi) above one revolution above required range
			//   if more than one revolution below required range, put it into the range that's [0, 2pi) below one revolution below required range
			// )
			if ( *this > upper_plus_one_rev ) {
				get_angle() = upper_plus_one_rev_in_rads + std::fmod( get_angle() - upper_plus_one_rev_in_rads, static_cast<T>( 2 * boost::math::constants::pi<double>() ) );
			}
			if ( *this < lower_less_one_rev ) {
				get_angle() = lower_less_one_rev_in_rads + std::fmod( get_angle() - lower_less_one_rev_in_rads, static_cast<T>( 2 * boost::math::constants::pi<double>() ) );
			}

			// Perform any necessary fine-tuning
			// This is more careful than the above code about how it handles the end-points
			while ( *this > upper || ( arg_endpoint_loc == angle_endpoint_loc::USE_LOWER && common::difference( *this, upper ) < endpoint_tolerance ) ) {
				( *this ) -= one_revolution<T>();
			}
			while ( *this < lower || ( arg_endpoint_loc == angle_endpoint_loc::USE_UPPER && common::difference( *this, lower ) < endpoint_tolerance ) ) {
				( *this ) += one_revolution<T>();
			}
			return *this;
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		angle<T> operator-(const angle<T> &arg_angle ///< TODOCUMENT
		                   ) {
			return ( zero_angle<T>() - arg_angle );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		T operator/(const angle<T> &arg_angle_a, ///< TODOCUMENT
		            const angle<T> &arg_angle_b  ///< TODOCUMENT
		            ) {
			return ( arg_angle_a.get_angle() / arg_angle_b.get_angle() );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		bool operator<(const angle<T> &arg_angle_a, ///< TODOCUMENT
		               const angle<T> &arg_angle_b  ///< TODOCUMENT
		               ) {
			return ( arg_angle_a.get_angle() < arg_angle_b.get_angle() );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		bool operator==(const angle<T> &arg_angle_a, ///< TODOCUMENT
		                const angle<T> &arg_angle_b  ///< TODOCUMENT
		                ) {
			return ( arg_angle_a.get_angle() == arg_angle_b.get_angle() );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		std::ostream & operator<<(std::ostream   &arg_os,   ///< TODOCUMENT
		                          const angle<T> &arg_angle ///< TODOCUMENT
		                          ) {
			// Splurge the whole thing out in one insertion operations
			// ( so the output can be positioned nicely with setw() )
			// (lexical_cast<string() was tried but it can result in excess precision)
			arg_os << ( std::to_string( angle_in_degrees( arg_angle ) ) + "°" ); // The degree symbol is ASCII character 176 (hex B0)
			return arg_os;
		}

		/// \brief Return the smallest difference between two angles when allowing for wrapping
		///        whilst assuming that the two angles are already shifted to the same range
		///
		/// This does less work than wrapped_difference() by assuming that both angles have already
		/// been shifted to the same range. If in doubt, use wrapped_difference().
		///
		/// \relates angle
		template <typename T>
		inline angle<T> unshifted_wrapped_difference(const angle<T> &arg_angle_a, ///< The first angle to compare
		                                             const angle<T> &arg_angle_b  ///< The second angle to compare
		                                             ) {
			const auto min_angle              = std::min( arg_angle_a, arg_angle_b );
			const auto max_angle              = std::max( arg_angle_a, arg_angle_b );
			const auto min_angle_plus_one_rev = min_angle + one_revolution<T>();
			return std::min(
				common::difference( arg_angle_a,            arg_angle_b ),
				common::difference( min_angle_plus_one_rev, max_angle   )
			);
		}

		/// \brief Return the smallest difference between two angles when allowing for wrapping
		///
		/// \relates angle
		template <typename T>
		inline angle<T> wrapped_difference(angle<T> arg_angle_a, ///< The first angle to compare
		                                   angle<T> arg_angle_b  ///< The second angle to compare
		                                   ) {
			arg_angle_a.shift();
			arg_angle_b.shift();
			return unshifted_wrapped_difference( arg_angle_a, arg_angle_b );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		angle<T> shift_copy(angle<T>                  arg_angle,                ///< TODOCUMENT
		                    const angle<T>           &arg_lower_angle_of_range, ///< TODOCUMENT
		                    const angle_endpoint_loc &arg_endpoint_loc          ///< TODOCUMENT
		                    ) {
			arg_angle.shift( arg_lower_angle_of_range, arg_endpoint_loc );
			return arg_angle;
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		inline angle<T> make_angle_from_degrees(const double &arg_degrees ///< TODOCUMENT
		                                        ) {
		//	return angle( quantity<plane_angle>( arg_degrees * degree ) );
			return angle<T>{ boost::numeric_cast<T>(
				( boost::math::constants::pi<double>() / 180.0 ) * arg_degrees
			) };
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		inline angle<T> make_angle_from_radians(const double &arg_radians ///< TODOCUMENT
		                                        ) {
		//	return angle( quantity<plane_angle>( arg_radians * radian ) );
			return angle<T>{ boost::numeric_cast<T>( arg_radians ) };
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		inline angle<T> make_angle_from_revolutions(const double &arg_revs ///< TODOCUMENT
		                                            ) {
		//	return angle( quantity<plane_angle>( arg_revs * revolution ) );
			return angle<T>{ boost::numeric_cast<T>( 2.0 * boost::math::constants::pi<double>() * arg_revs ) };
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		T angle_in_degrees(const angle<T> &arg_angle ///< TODOCUMENT
		                   ) {
		//	return quantity<units::degree::plane_angle>( arg_angle.get_angle() ).value();
			return static_cast<T>( 180.0 / boost::math::constants::pi<double>() ) * arg_angle.get_angle();
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		T angle_in_radians(const angle<T> &arg_angle ///< TODOCUMENT
		                   ) {
		//	return arg_angle.get_angle().value();
			return arg_angle.get_angle();
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		T angle_in_revolutions(const angle<T> &arg_angle ///< TODOCUMENT
		                       ) {
		//	return quantity<units::revolution::plane_angle>( arg_angle.get_angle() ).value();
			return ( arg_angle.get_angle() / ( 2.0 * boost::math::constants::pi<double>() ) );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		angle<T> zero_angle() {
			return make_angle_from_radians<T>( 0.0 );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		angle<T> half_revolution() {
			return make_angle_from_revolutions<T>( 0.5 );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		angle<T> one_revolution() {
			return make_angle_from_revolutions<T>( 1.0 );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T, typename U>
		angle<T> convert_angle_type(const angle<U> &arg_angle
		                            ) {
			return make_angle_from_radians<T>( angle_in_radians( arg_angle ) );
		}

		/// \brief TODOCUMENT
		using float_angle = angle<float>;

		/// \brief TODOCUMENT
		using doub_angle  = angle<double>;
	} // namespace geom
} // namespace cath

#endif
