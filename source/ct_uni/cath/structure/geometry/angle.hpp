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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ANGLE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ANGLE_HPP

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

namespace cath::geom {
	template <typename T>
	class angle;

	template<typename T>
	constexpr angle<T> make_angle_from_degrees(const double &);


	/// \brief TODOCUMENT
	template <typename T>
	class angle final : private boost::totally_ordered<angle<T>,
	                            boost::multiplicative <angle<T>, T> > {
	private:
		using angle_type = T;
//		using angle_type = boost::units::quantity<boost::units::si::plane_angle>;

		/// \brief TODOCUMENT
		angle_type angle_value;

		constexpr angle_type & get_angle();

	public:
		explicit constexpr angle(const angle_type &);

		constexpr const angle_type & get_angle() const;

		constexpr angle & operator-=(const angle &);
		constexpr angle & operator+=(const angle &);
		constexpr angle & operator*=(const T &);
		constexpr angle & operator/=(const T &);

		constexpr angle & quick_shift();

		constexpr angle & shift(const angle & = make_angle_from_degrees<angle_type>( 0.0 ),
		                        const angle_endpoint_loc & = angle_endpoint_loc::USE_LOWER);
	};

	/// \brief TODOCUMENT
	///
	/// \relates angle
	///
	/// \param prm_degrees TODOCUMENT
	template <typename T>
	constexpr angle<T> make_angle_from_degrees( const double &prm_degrees ) {
		return angle<T>{ static_cast<T>( ( boost::math::constants::pi<double>() / 180.0 ) * prm_degrees ) };
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	///
	/// \param prm_radians TODOCUMENT
	template <typename T>
	constexpr angle<T> make_angle_from_radians( const double &prm_radians ) {
		return angle<T>{ static_cast<T>( prm_radians ) };
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	///
	/// \param prm_revs TODOCUMENT
	template <typename T>
	constexpr angle<T> make_angle_from_revolutions( const double &prm_revs ) {
		return angle<T>{ static_cast<T>( 2.0 * boost::math::constants::pi<double>() * prm_revs ) };
	}

	/// \brief TODOCUMENT
	template <typename T>
	inline constexpr angle<T> ZERO_ANGLE = make_angle_from_radians<T>( 0.0 );

	/// \brief TODOCUMENT
	template <typename T>
	inline constexpr angle<T> HALF_REVOLUTION = make_angle_from_revolutions<T>( 0.5 );

	/// \brief TODOCUMENT
	template <typename T>
	inline constexpr angle<T> ONE_REVOLUTION = make_angle_from_revolutions<T>( 1.0 );


	/// \brief Ctor for angle
	///
	/// \param prm_angle_value TODOCUMENT
	template <typename T>
	constexpr angle<T>::angle( const angle_type &prm_angle_value ) : angle_value( prm_angle_value ) {
		if ( ! ( prm_angle_value > ::std::numeric_limits<angle_type>::lowest()
		     && prm_angle_value < ::std::numeric_limits<angle_type>::max() ) ) {
			BOOST_THROW_EXCEPTION( cath::common::invalid_argument_exception( "Angle value must be finite" ) );
		}
	}

	/// \brief TODOCUMENT
	template <typename T>
	constexpr typename angle<T>::angle_type & angle<T>::get_angle() {
		return angle_value;
	}

	/// \brief TODOCUMENT
	template <typename T>
	constexpr const typename angle<T>::angle_type & angle<T>::get_angle() const {
		return angle_value;
	}

	/// \brief TODOCUMENT
	template <typename T>
	constexpr angle<T> & angle<T>::operator-=(const angle<T> &prm_angle ///< TODOCUMENT
	                                       ) {
		get_angle() -= prm_angle.get_angle();
		return *this;
	}

	/// \brief TODOCUMENT
	template <typename T>
	constexpr angle<T> & angle<T>::operator+=(const angle<T> &prm_angle ///< TODOCUMENT
	                                          ) {
		get_angle() += prm_angle.get_angle();
		return *this;
	}

	/// \brief TODOCUMENT
	template <typename T>
	constexpr angle<T> & angle<T>::operator*=(const T &prm_angle ///< TODOCUMENT
	                                          ) {
		get_angle() *= prm_angle;
		return *this;
	}

	/// \brief TODOCUMENT
	template <typename T>
	constexpr angle<T> & angle<T>::operator/=(const T &prm_angle ///< TODOCUMENT
	                                          ) {
		get_angle() /= prm_angle;
		return *this;
	}

	/// \brief TODOCUMENT
	template <typename T>
	constexpr angle<T> & angle<T>::quick_shift() {
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
	constexpr angle<T> & angle<T>::shift(const angle<T>           &prm_lower_angle_of_range, ///< TODOCUMENT
	                                     const angle_endpoint_loc &prm_endpoint_loc          ///< TODOCUMENT
	                                     ) {
		// Grab some basics
		const angle  endpoint_tolerance         = make_angle_from_degrees<T>( static_cast<T>( 0.0000000001 ) );
		const angle &lower                      = prm_lower_angle_of_range;
		const angle  upper                      = prm_lower_angle_of_range + ONE_REVOLUTION<T>;
		const angle  lower_less_one_rev         = lower - ONE_REVOLUTION<T>;
		const angle  upper_plus_one_rev         = upper + ONE_REVOLUTION<T>;
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
			( *this ) -= ONE_REVOLUTION<T>;
		}
		while ( *this < lower || ( prm_endpoint_loc == angle_endpoint_loc::USE_UPPER && common::difference( *this, lower ) < endpoint_tolerance ) ) {
			( *this ) += ONE_REVOLUTION<T>;
		}
		return *this;
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	///
	/// \param prm_lhs TODOCUMENT
	/// \param prm_Rhs TODOCUMENT
	template <typename T>
	constexpr angle<T> operator+( angle<T> prm_lhs, const angle<T> &prm_rhs ) {
		prm_lhs += prm_rhs;
		return prm_lhs;
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	///
	/// \param prm_lhs TODOCUMENT
	/// \param prm_Rhs TODOCUMENT
	template <typename T>
	constexpr angle<T> operator-( angle<T> prm_lhs, const angle<T> &prm_rhs ) {
		prm_lhs -= prm_rhs;
		return prm_lhs;
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	template <typename T>
	constexpr angle<T> operator-(const angle<T> &prm_angle ///< TODOCUMENT
	                             ) {
		return ( ZERO_ANGLE<T> - prm_angle );
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	template <typename T>
	constexpr T operator/(const angle<T> &prm_angle_a, ///< TODOCUMENT
	                      const angle<T> &prm_angle_b  ///< TODOCUMENT
	                      ) {
		return ( prm_angle_a.get_angle() / prm_angle_b.get_angle() );
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	template <typename T>
	constexpr bool operator<(const angle<T> &prm_angle_a, ///< TODOCUMENT
	                         const angle<T> &prm_angle_b  ///< TODOCUMENT
	                         ) {
		return ( prm_angle_a.get_angle() < prm_angle_b.get_angle() );
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	template <typename T>
	constexpr bool operator==(const angle<T> &prm_angle_a, ///< TODOCUMENT
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
	constexpr angle<T> unshifted_wrapped_difference(const angle<T> &prm_angle_a, ///< The first angle to compare
	                                                const angle<T> &prm_angle_b  ///< The second angle to compare
	                                                ) {
		const auto min_angle              = std::min( prm_angle_a, prm_angle_b );
		const auto max_angle              = std::max( prm_angle_a, prm_angle_b );
		const auto min_angle_plus_one_rev = min_angle + ONE_REVOLUTION<T>;
		return std::min(
			common::difference( prm_angle_a,            prm_angle_b ),
			common::difference( min_angle_plus_one_rev, max_angle   )
		);
	}

	/// \brief Return the smallest difference between two angles when allowing for wrapping
	///
	/// \relates angle
	template <typename T>
	constexpr angle<T> wrapped_difference(angle<T> prm_angle_a, ///< The first angle to compare
	                                      angle<T> prm_angle_b  ///< The second angle to compare
	                                      ) {
		prm_angle_a.shift();
		prm_angle_b.shift();
		return unshifted_wrapped_difference( prm_angle_a, prm_angle_b );
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	///
	/// \param prm_angle                TODOCUMENT
	/// \param prm_lower_angle_of_range TODOCUMENT
	/// \param prm_endpoint_loc         TODOCUMENT
	template <typename T>
	constexpr angle<T> shift_copy( angle<T>        prm_angle,
	                               const angle<T> &prm_lower_angle_of_range   = make_angle_from_degrees<T>( 0.0 ),
	                               const angle_endpoint_loc &prm_endpoint_loc = angle_endpoint_loc::USE_LOWER ) {
		prm_angle.shift( prm_lower_angle_of_range, prm_endpoint_loc );
		return prm_angle;
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	template <typename T>
	constexpr T angle_in_degrees(const angle<T> &prm_angle ///< TODOCUMENT
	                             ) {
	//	return quantity<units::degree::plane_angle>( prm_angle.get_angle() ).value();
		return static_cast<T>( 180.0 / boost::math::constants::pi<double>() ) * prm_angle.get_angle();
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	template <typename T>
	constexpr T angle_in_radians(const angle<T> &prm_angle ///< TODOCUMENT
	                             ) {
	//	return prm_angle.get_angle().value();
		return prm_angle.get_angle();
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	template <typename T>
	constexpr T angle_in_revolutions(const angle<T> &prm_angle ///< TODOCUMENT
	                                 ) {
	//	return quantity<units::revolution::plane_angle>( prm_angle.get_angle() ).value();
		return ( prm_angle.get_angle() / ( 2.0 * boost::math::constants::pi<double>() ) );
	}

	/// \brief TODOCUMENT
	///
	/// \relates angle
	template <typename T, typename U>
	constexpr angle<T> convert_angle_type(const angle<U> &prm_angle
	                            ) {
		return make_angle_from_radians<T>( angle_in_radians( prm_angle ) );
	}

	/// \brief TODOCUMENT
	using float_angle = angle<float>;

	/// \brief TODOCUMENT
	using doub_angle  = angle<double>;

} // namespace cath::geom

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ANGLE_HPP
