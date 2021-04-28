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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ANGLE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ANGLE_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/operators.hpp>
//#include <boost/units/quantity.hpp>
//#include <boost/units/systems/si/plane_angle.hpp>

#include "cath/common/difference.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/structure/geometry/angle_endpoint_loc.hpp"

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
		inline angle<T>::angle(const angle_type &prm_angle_value ///< TODOCUMENT
		                       ) : angle_value( prm_angle_value ) {
			if ( ! ::boost::math::isfinite( prm_angle_value ) ) {
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
		inline angle<T> & angle<T>::operator-=(const angle<T> &prm_angle ///< TODOCUMENT
		                                       ) {
			get_angle() -= prm_angle.get_angle();
			return *this;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline angle<T> & angle<T>::operator+=(const angle<T> &prm_angle ///< TODOCUMENT
		                                       ) {
			get_angle() += prm_angle.get_angle();
			return *this;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline angle<T> & angle<T>::operator*=(const T &prm_angle ///< TODOCUMENT
		                                       ) {
			get_angle() *= prm_angle;
			return *this;
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline angle<T> & angle<T>::operator/=(const T &prm_angle ///< TODOCUMENT
		                                       ) {
			get_angle() /= prm_angle;
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
		inline angle<T> & angle<T>::shift(const angle<T>           &prm_lower_angle_of_range, ///< TODOCUMENT
		                                  const angle_endpoint_loc &prm_endpoint_loc          ///< TODOCUMENT
		                                  ) {
			// Grab some basics
			const angle  endpoint_tolerance         = make_angle_from_degrees<T>( static_cast<T>( 0.0000000001 ) );
			const angle &lower                      = prm_lower_angle_of_range;
			const angle  upper                      = prm_lower_angle_of_range + one_revolution<T>();
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
			while ( *this > upper || ( prm_endpoint_loc == angle_endpoint_loc::USE_LOWER && common::difference( *this, upper ) < endpoint_tolerance ) ) {
				( *this ) -= one_revolution<T>();
			}
			while ( *this < lower || ( prm_endpoint_loc == angle_endpoint_loc::USE_UPPER && common::difference( *this, lower ) < endpoint_tolerance ) ) {
				( *this ) += one_revolution<T>();
			}
			return *this;
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		angle<T> operator-(const angle<T> &prm_angle ///< TODOCUMENT
		                   ) {
			return ( zero_angle<T>() - prm_angle );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		T operator/(const angle<T> &prm_angle_a, ///< TODOCUMENT
		            const angle<T> &prm_angle_b  ///< TODOCUMENT
		            ) {
			return ( prm_angle_a.get_angle() / prm_angle_b.get_angle() );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		bool operator<(const angle<T> &prm_angle_a, ///< TODOCUMENT
		               const angle<T> &prm_angle_b  ///< TODOCUMENT
		               ) {
			return ( prm_angle_a.get_angle() < prm_angle_b.get_angle() );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		bool operator==(const angle<T> &prm_angle_a, ///< TODOCUMENT
		                const angle<T> &prm_angle_b  ///< TODOCUMENT
		                ) {
			return ( prm_angle_a.get_angle() == prm_angle_b.get_angle() );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		std::ostream & operator<<(std::ostream   &prm_os,   ///< TODOCUMENT
		                          const angle<T> &prm_angle ///< TODOCUMENT
		                          ) {
			// Splurge the whole thing out in one insertion operations
			// ( so the output can be positioned nicely with setw() )
			// (lexical_cast<string() was tried but it can result in excess precision)
			prm_os << ( std::to_string( angle_in_degrees( prm_angle ) ) + "°" ); // The degree symbol is ASCII character 176 (hex B0)
			return prm_os;
		}

		/// \brief Return the smallest difference between two angles when allowing for wrapping
		///        whilst assuming that the two angles are already shifted to the same range
		///
		/// This does less work than wrapped_difference() by assuming that both angles have already
		/// been shifted to the same range. If in doubt, use wrapped_difference().
		///
		/// \relates angle
		template <typename T>
		inline angle<T> unshifted_wrapped_difference(const angle<T> &prm_angle_a, ///< The first angle to compare
		                                             const angle<T> &prm_angle_b  ///< The second angle to compare
		                                             ) {
			const auto min_angle              = std::min( prm_angle_a, prm_angle_b );
			const auto max_angle              = std::max( prm_angle_a, prm_angle_b );
			const auto min_angle_plus_one_rev = min_angle + one_revolution<T>();
			return std::min(
				common::difference( prm_angle_a,            prm_angle_b ),
				common::difference( min_angle_plus_one_rev, max_angle   )
			);
		}

		/// \brief Return the smallest difference between two angles when allowing for wrapping
		///
		/// \relates angle
		template <typename T>
		inline angle<T> wrapped_difference(angle<T> prm_angle_a, ///< The first angle to compare
		                                   angle<T> prm_angle_b  ///< The second angle to compare
		                                   ) {
			prm_angle_a.shift();
			prm_angle_b.shift();
			return unshifted_wrapped_difference( prm_angle_a, prm_angle_b );
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		angle<T> shift_copy(angle<T>                  prm_angle,                ///< TODOCUMENT
		                    const angle<T>           &prm_lower_angle_of_range, ///< TODOCUMENT
		                    const angle_endpoint_loc &prm_endpoint_loc          ///< TODOCUMENT
		                    ) {
			prm_angle.shift( prm_lower_angle_of_range, prm_endpoint_loc );
			return prm_angle;
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		inline angle<T> make_angle_from_degrees(const double &prm_degrees ///< TODOCUMENT
		                                        ) {
		//	return angle( quantity<plane_angle>( prm_degrees * degree ) );
			return angle<T>{ boost::numeric_cast<T>(
				( boost::math::constants::pi<double>() / 180.0 ) * prm_degrees
			) };
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		inline angle<T> make_angle_from_radians(const double &prm_radians ///< TODOCUMENT
		                                        ) {
		//	return angle( quantity<plane_angle>( prm_radians * radian ) );
			return angle<T>{ boost::numeric_cast<T>( prm_radians ) };
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		inline angle<T> make_angle_from_revolutions(const double &prm_revs ///< TODOCUMENT
		                                            ) {
		//	return angle( quantity<plane_angle>( prm_revs * revolution ) );
			return angle<T>{ boost::numeric_cast<T>( 2.0 * boost::math::constants::pi<double>() * prm_revs ) };
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		T angle_in_degrees(const angle<T> &prm_angle ///< TODOCUMENT
		                   ) {
		//	return quantity<units::degree::plane_angle>( prm_angle.get_angle() ).value();
			return static_cast<T>( 180.0 / boost::math::constants::pi<double>() ) * prm_angle.get_angle();
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		T angle_in_radians(const angle<T> &prm_angle ///< TODOCUMENT
		                   ) {
		//	return prm_angle.get_angle().value();
			return prm_angle.get_angle();
		}

		/// \brief TODOCUMENT
		///
		/// \relates angle
		template <typename T>
		T angle_in_revolutions(const angle<T> &prm_angle ///< TODOCUMENT
		                       ) {
		//	return quantity<units::revolution::plane_angle>( prm_angle.get_angle() ).value();
			return ( prm_angle.get_angle() / ( 2.0 * boost::math::constants::pi<double>() ) );
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
		angle<T> convert_angle_type(const angle<U> &prm_angle
		                            ) {
			return make_angle_from_radians<T>( angle_in_radians( prm_angle ) );
		}

		/// \brief TODOCUMENT
		using float_angle = angle<float>;

		/// \brief TODOCUMENT
		using doub_angle  = angle<double>;
	} // namespace geom
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ANGLE_HPP
