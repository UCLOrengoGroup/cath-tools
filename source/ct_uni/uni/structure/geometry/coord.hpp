/// \file
/// \brief The coord class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_GEOMETRY_COORD_HPP
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_GEOMETRY_COORD_HPP

#include <boost/geometry/geometries/geometries.hpp> // ***** TEMPORARY? *****
#include <boost/geometry/geometries/point.hpp> // ***** TEMPORARY? *****
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/operators.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

#include "common/detail/maybe_unused_namespace_scope_constexpr.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/json_style.hpp"
#include "common/property_tree/read_from_ptree.hpp"

namespace cath { namespace geom { template <typename T> class angle; } }
namespace cath { namespace geom { using doub_angle = angle<double>; } }

namespace cath {
	namespace geom {

		/// \brief TODOCUMENT
		///
		/// \todo Seriously consider implementing these geometry-based classes using Boost geometry:
		///       www.boost.org/libs/geometry/
		/// \todo Consider making each axis into an enum value and changing the getters to one that
		///       requires an axis argument (and make get_x(),get_y() and get_z() into
		///       non-member, non-friend functions that call that getter).
		class coord final : boost::equality_comparable<coord> {
		private:
			/// \brief TODOCUMENT
			double x;

			/// \brief TODOCUMENT
			double y;

			/// \brief TODOCUMENT
			double z;

		public:
			/// \brief TODOCUMENT
			static constexpr size_t NUM_DIMS = 3;

			/// \brief TODOCUMENT
			static constexpr double TOLERANCE_FOR_COORD_CLOSENESS_CHECKS = 0.00001;

			constexpr coord(const double &,
			                const double &,
			                const double &);

			template <typename T>
			explicit coord(const boost::geometry::model::point<T, 3, boost::geometry::cs::cartesian> &);

			constexpr const double & get_x() const;
			constexpr const double & get_y() const;
			constexpr const double & get_z() const;

			coord & operator-=(const coord &);
			coord & operator+=(const coord &);
			coord & operator*=(const double &);
			coord & operator/=(const double &);

			/// \brief TODOCUMENT
			///
			/// \todo Come Boost 1.58.0 (or C++17, if Herb Sutter has gotten his way (n4029)), just use braced list here
			template <typename T>
			operator boost::geometry::model::point<T, 3, boost::geometry::cs::cartesian>() const {
				return boost::geometry::model::point<T, 3, boost::geometry::cs::cartesian>{
					boost::numeric_cast<T>( get_x() ),
					boost::numeric_cast<T>( get_y() ),
					boost::numeric_cast<T>( get_z() )
				};
			}
		};

		/// \brief Non-default constructor for coord.
		inline constexpr coord::coord(const double &prm_x, ///< TODOCUMENT
		                              const double &prm_y, ///< TODOCUMENT
		                              const double &prm_z  ///< TODOCUMENT
		                              ) : x ( prm_x ),
		                                  y ( prm_y ),
		                                  z ( prm_z ) {
	// #ifndef NDEBUG
	// 		if ( ! boost::math::isfinite( prm_x ) || ! boost::math::isfinite( prm_y ) || ! boost::math::isfinite( prm_z ) ) {
	// 			BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Arguments x, y and z must be a normal, finite floating-point numbers"));
	// 		}
	// #endif
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline coord::coord(const boost::geometry::model::point<T, 3, boost::geometry::cs::cartesian> &prm_point /// \brief TODOCUMENT
		                    ) : coord(
		                        	boost::numeric_cast<double>( prm_point.template get<0>() ),
		                        	boost::numeric_cast<double>( prm_point.template get<1>() ),
		                        	boost::numeric_cast<double>( prm_point.template get<2>() )
		                        ) {
		}

		/// \brief TODOCUMENT
		inline constexpr const double & coord::get_x() const {
			return x;
		}

		/// \brief TODOCUMENT
		inline constexpr const double & coord::get_y() const {
			return y;
		}

		/// \brief TODOCUMENT
		inline constexpr const double & coord::get_z() const {
			return z;
		}

		/// \brief TODOCUMENT
		inline coord & coord::operator-=(const coord &prm_coord ///< TODOCUMENT
		                                 ) {
			x -= prm_coord.get_x();
			y -= prm_coord.get_y();
			z -= prm_coord.get_z();
			return ( *this );
		}

		/// \brief TODOCUMENT
		inline coord & coord::operator+=(const coord &prm_coord ///< TODOCUMENT
		                                  ) {
			x += prm_coord.get_x();
			y += prm_coord.get_y();
			z += prm_coord.get_z();
			return ( *this );
		}

		/// \brief TODOCUMENT
		inline coord & coord::operator*=(const double &prm_factor ///< TODOCUMENT
		                                 ) {
	#ifndef NDEBUG
			if ( boost::math::isnan( prm_factor ) ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("InvalidcoordFactor"));
			}
	#endif
			x *= prm_factor;
			y *= prm_factor;
			z *= prm_factor;
			return ( *this );
		}

		/// \brief TODOCUMENT
		inline coord & coord::operator/=(const double &prm_factor ///< TODOCUMENT
		                                 ) {
	#ifndef NDEBUG
			if ( ! boost::math::isnormal(prm_factor) ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("InvalidcoordFactor"));
			}
	#endif
			return operator*=( 1.0 / prm_factor );
		}

		/// \brief TODOCUMENT
		constexpr coord  ORIGIN_COORD { 0.0, 0.0, 0.0 };

		/// \brief TODOCUMENT
		constexpr coord  UNIT_X_COORD { 1.0, 0.0, 0.0 };
		MAYBE_UNUSED_NAMESPACE_SCOPE_CONSTEXPR( UNIT_X_COORD )

		/// \brief TODOCUMENT
		constexpr coord  UNIT_Y_COORD { 0.0, 1.0, 0.0 };
		MAYBE_UNUSED_NAMESPACE_SCOPE_CONSTEXPR( UNIT_Y_COORD )

		/// \brief TODOCUMENT
		constexpr coord  UNIT_Z_COORD { 0.0, 0.0, 1.0 };
		MAYBE_UNUSED_NAMESPACE_SCOPE_CONSTEXPR( UNIT_Z_COORD )


		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::additive<>
		///
		/// \relates coord
		inline constexpr coord operator-(const coord &prm_coord_lhs, ///< TODOCUMENT
		                                 const coord &prm_coord_rhs  ///< TODOCUMENT
		                                 ) {
			return {
				prm_coord_lhs.get_x() - prm_coord_rhs.get_x(),
				prm_coord_lhs.get_y() - prm_coord_rhs.get_y(),
				prm_coord_lhs.get_z() - prm_coord_rhs.get_z(),
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::additive<>
		///
		/// \relates coord
		inline constexpr coord operator+(const coord &prm_coord_lhs, ///< TODOCUMENT
		                                 const coord &prm_coord_rhs  ///< TODOCUMENT
		                                 ) {
			return {
				prm_coord_lhs.get_x() + prm_coord_rhs.get_x(),
				prm_coord_lhs.get_y() + prm_coord_rhs.get_y(),
				prm_coord_lhs.get_z() + prm_coord_rhs.get_z(),
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::multiplicative<>
		///
		/// \relates coord
		inline constexpr coord operator*(const coord  &prm_coord, ///< TODOCUMENT
		                                 const double &prm_factor ///< TODOCUMENT
		                                 ) {
			return {
				prm_coord.get_x() * prm_factor,
				prm_coord.get_y() * prm_factor,
				prm_coord.get_z() * prm_factor
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::multiplicative<>
		///
		/// \relates coord
		inline constexpr coord operator*(const double &prm_factor, ///< TODOCUMENT
		                                 const coord  &prm_coord   ///< TODOCUMENT
		                                 ) {
			return {
				prm_coord.get_x() * prm_factor,
				prm_coord.get_y() * prm_factor,
				prm_coord.get_z() * prm_factor
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::multiplicative<>
		///
		/// \relates coord
		inline constexpr coord operator/(const coord  &prm_coord, ///< TODOCUMENT
		                                 const double &prm_factor ///< TODOCUMENT
		                                 ) {
			return {
				prm_coord.get_x() / prm_factor,
				prm_coord.get_y() / prm_factor,
				prm_coord.get_z() / prm_factor
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::multiplicative<>
		///
		/// \relates coord
		inline constexpr coord operator/(const double &prm_factor, ///< TODOCUMENT
		                                 const coord  &prm_coord   ///< TODOCUMENT
		                                 ) {
			return {
				prm_coord.get_x() / prm_factor,
				prm_coord.get_y() / prm_factor,
				prm_coord.get_z() / prm_factor
			};
		}




		constexpr bool operator==(const coord &,
		                          const coord &);
		coord operator-(const coord &);

		bool not_zero(const coord &);
		double dot_product(const coord &,
		                   const coord &);
		coord cross_product(const coord &,
		                    const coord &);
		coord int_cast_copy(const coord &);

		coord normalise_copy(const coord &);
		coord parallel_component_copy(const coord &,
		                              const coord &);
		coord perpendicular_component_copy(const coord &,
		                                   const coord &);
		double length(const coord &);
		double distance_between_points(const coord &,
		                               const coord &);

		doub_angle angle_between_two_vectors(const coord &,
		                                     const coord &);

		doub_angle angle_between_three_points(const coord &,
		                                      const coord &,
		                                      const coord &);

		doub_angle dihedral_angle_between_four_points(const coord &,
		                                              const coord &,
		                                              const coord &,
		                                              const coord &);

		template <typename T>
		angle<T> angle_between_two_vectors(const coord &,
		                                   const coord &);

		template <typename T>
		angle<T> angle_between_three_points(const coord &,
		                                    const coord &,
		                                    const coord &);

		template <typename T>
		angle<T> dihedral_angle_between_four_points(const coord &,
		                                            const coord &,
		                                            const coord &,
		                                            const coord &);

		doub_angle planar_angle_between(const coord &,
		                                const coord &,
		                                const coord &);

		std::string to_string(const coord &);
		std::ostream & operator<<(std::ostream &,
		                          const coord &);

		coord coord_from_ptree(const boost::property_tree::ptree &);
		void save_to_ptree(boost::property_tree::ptree &,
		                   const coord &);

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline constexpr double squared_length(const coord &prm_coord ///< TODOCUMENT
		                                       ) {
			return (
				prm_coord.get_x() * prm_coord.get_x() +
				prm_coord.get_y() * prm_coord.get_y() +
				prm_coord.get_z() * prm_coord.get_z()
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline double length(const coord &prm_coord ///< TODOCUMENT
		                     ) {
			return sqrt( squared_length( prm_coord ) );
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline double distance_between_points(const coord &prm_coord1, ///< TODOCUMENT
		                                      const coord &prm_coord2  ///< TODOCUMENT
		                                      ) {
			// Equivalently, either :
			//   length(prm_coord2 - prm_coord1)
			// or
			//   sqrt(squared_distance_between_points(prm_coord1, prm_coord2));
			return length( prm_coord2 - prm_coord1 );
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline constexpr double squared_distance_between_points(const coord &prm_coord1, ///< TODOCUMENT
		                                                        const coord &prm_coord2  ///< TODOCUMENT
		                                                        ) {
			return squared_length( prm_coord2 - prm_coord1 );
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline constexpr bool operator==(const coord &prm_coord1, ///< TODOCUMENT
		                                 const coord &prm_coord2  ///< TODOCUMENT
		                                 ) {
			return (
				squared_distance_between_points( prm_coord1, prm_coord2 )
				<
				coord::TOLERANCE_FOR_COORD_CLOSENESS_CHECKS * coord::TOLERANCE_FOR_COORD_CLOSENESS_CHECKS
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline coord operator-(const coord &prm_coord ///< TODOCUMENT
		                       ) {
			return ( ORIGIN_COORD - prm_coord );
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline bool not_zero(const coord &prm_coord ///< TODOCUMENT
		                     ) {
			if ( prm_coord.get_x() != 0.0 ) {
				return true;
			}
			if ( prm_coord.get_y() != 0.0 ) {
				return true;
			}
			if ( prm_coord.get_z() != 0.0 ) {
				return true;
			}
			return false;
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline double dot_product(const coord &prm_coord1, ///< TODOCUMENT
		                          const coord &prm_coord2  ///< TODOCUMENT
		                          ) {
			return (
				prm_coord1.get_x() * prm_coord2.get_x() +
				prm_coord1.get_y() * prm_coord2.get_y() +
				prm_coord1.get_z() * prm_coord2.get_z()
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline coord cross_product(const coord &prm_coord1, ///< TODOCUMENT
		                           const coord &prm_coord2  ///< TODOCUMENT
		                           ) {
			return coord(
				prm_coord1.get_y() * prm_coord2.get_z() - prm_coord1.get_z() * prm_coord2.get_y(),
				prm_coord1.get_z() * prm_coord2.get_x() - prm_coord1.get_x() * prm_coord2.get_z(),
				prm_coord1.get_x() * prm_coord2.get_y() - prm_coord1.get_y() * prm_coord2.get_x()
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline coord int_cast_copy(const coord &prm_coord ///< TODOCUMENT
		                           ) {
			return coord(
				boost::numeric_cast<int>( prm_coord.get_x() ),
				boost::numeric_cast<int>( prm_coord.get_y() ),
				boost::numeric_cast<int>( prm_coord.get_z() )
			);
			return prm_coord;
		}

		/// \brief Return the scalar component of the first specified coord in the direction of the second
		///
		/// \relates coord
		inline double scalar_component(const coord &prm_source_coord, ///< The coord of which to return the scalar component
		                               const coord &prm_dirn          ///< The direction in which the component should be taken (does not have to be normalised)
		                               ) {
			return dot_product( prm_source_coord, normalise_copy( prm_dirn ) );
		}

		/// \brief Return the parallel component of the first specified coord in the direction of the second
		///
		/// \relates coord
		inline coord parallel_component_copy(const coord &prm_source_coord, ///< The coord of which to return the parallel component
		                                     const coord &prm_dirn          ///< The direction in which the component should be taken (does not have to be normalised)
		                                     ) {
			return scalar_component( prm_source_coord, prm_dirn ) * normalise_copy( prm_dirn );
		}

		/// \brief Return the perpendicular component of the first specified coord in the direction perpendicular to the second
		///
		/// \relates coord
		inline coord perpendicular_component_copy(const coord &prm_source_coord, ///< The coord of which to return the perpendicular component
		                                          const coord &prm_dirn          ///< The direction perpendicular to which the component should be taken (does not have to be normalised)
		                                          ) {
			const coord  unit_dirn = normalise_copy( prm_dirn );
			const double factor    = dot_product   ( prm_source_coord, unit_dirn );
			return prm_source_coord - ( factor * unit_dirn );
		}


	} // namespace geom
	
	namespace common {
	
		/// \brief Specialisation of cath::common::read_from_ptree for coord
		template <>
		inline geom::coord read_from_ptree<geom::coord>(const boost::property_tree::ptree &prm_ptree ///< The ptree from which to read the coord
		                                                ) {
			return geom::coord_from_ptree( prm_ptree );
		}
	
	} // namespace common
} // namespace cath

#endif
