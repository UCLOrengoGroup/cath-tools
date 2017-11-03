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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_COORD_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_COORD_H

#include <boost/geometry/geometries/geometries.hpp> // ***** TEMPORARY? *****
#include <boost/geometry/geometries/point.hpp> // ***** TEMPORARY? *****
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/operators.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

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
			static const coord  ORIGIN_COORD;

			/// \brief TODOCUMENT
			static const coord  UNIT_X;

			/// \brief TODOCUMENT
			static const coord  UNIT_Y;

			/// \brief TODOCUMENT
			static const coord  UNIT_Z;

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
		inline constexpr coord::coord(const double &arg_x, ///< TODOCUMENT
		                              const double &arg_y, ///< TODOCUMENT
		                              const double &arg_z  ///< TODOCUMENT
		                              ) : x ( arg_x ),
		                                  y ( arg_y ),
		                                  z ( arg_z ) {
	// #ifndef NDEBUG
	// 		if ( ! boost::math::isfinite( arg_x ) || ! boost::math::isfinite( arg_y ) || ! boost::math::isfinite( arg_z ) ) {
	// 			BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Arguments x, y and z must be a normal, finite floating-point numbers"));
	// 		}
	// #endif
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline coord::coord(const boost::geometry::model::point<T, 3, boost::geometry::cs::cartesian> &arg_point /// \brief TODOCUMENT
		                    ) : coord(
		                        	boost::numeric_cast<double>( arg_point.template get<0>() ),
		                        	boost::numeric_cast<double>( arg_point.template get<1>() ),
		                        	boost::numeric_cast<double>( arg_point.template get<2>() )
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
		inline coord & coord::operator-=(const coord &arg_coord ///< TODOCUMENT
		                                 ) {
			x -= arg_coord.get_x();
			y -= arg_coord.get_y();
			z -= arg_coord.get_z();
			return ( *this );
		}

		/// \brief TODOCUMENT
		inline coord & coord::operator+=(const coord &arg_coord ///< TODOCUMENT
		                                  ) {
			x += arg_coord.get_x();
			y += arg_coord.get_y();
			z += arg_coord.get_z();
			return ( *this );
		}

		/// \brief TODOCUMENT
		inline coord & coord::operator*=(const double &arg_factor ///< TODOCUMENT
		                                 ) {
	#ifndef NDEBUG
			if ( boost::math::isnan( arg_factor ) ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("InvalidcoordFactor"));
			}
	#endif
			x *= arg_factor;
			y *= arg_factor;
			z *= arg_factor;
			return ( *this );
		}

		/// \brief TODOCUMENT
		inline coord & coord::operator/=(const double &arg_factor ///< TODOCUMENT
		                                 ) {
	#ifndef NDEBUG
			if ( ! boost::math::isnormal(arg_factor) ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("InvalidcoordFactor"));
			}
	#endif
			return operator*=( 1.0 / arg_factor );
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::additive<>
		///
		/// \relates coord
		inline constexpr coord operator-(const coord &arg_coord_lhs, ///< TODOCUMENT
		                                 const coord &arg_coord_rhs  ///< TODOCUMENT
		                                 ) {
			return {
				arg_coord_lhs.get_x() - arg_coord_rhs.get_x(),
				arg_coord_lhs.get_y() - arg_coord_rhs.get_y(),
				arg_coord_lhs.get_z() - arg_coord_rhs.get_z(),
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::additive<>
		///
		/// \relates coord
		inline constexpr coord operator+(const coord &arg_coord_lhs, ///< TODOCUMENT
		                                 const coord &arg_coord_rhs  ///< TODOCUMENT
		                                 ) {
			return {
				arg_coord_lhs.get_x() + arg_coord_rhs.get_x(),
				arg_coord_lhs.get_y() + arg_coord_rhs.get_y(),
				arg_coord_lhs.get_z() + arg_coord_rhs.get_z(),
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::multiplicative<>
		///
		/// \relates coord
		inline constexpr coord operator*(const coord  &arg_coord, ///< TODOCUMENT
		                                 const double &arg_factor ///< TODOCUMENT
		                                 ) {
			return {
				arg_coord.get_x() * arg_factor,
				arg_coord.get_y() * arg_factor,
				arg_coord.get_z() * arg_factor
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::multiplicative<>
		///
		/// \relates coord
		inline constexpr coord operator*(const double &arg_factor, ///< TODOCUMENT
		                                 const coord  &arg_coord   ///< TODOCUMENT
		                                 ) {
			return {
				arg_coord.get_x() * arg_factor,
				arg_coord.get_y() * arg_factor,
				arg_coord.get_z() * arg_factor
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::multiplicative<>
		///
		/// \relates coord
		inline constexpr coord operator/(const coord  &arg_coord, ///< TODOCUMENT
		                                 const double &arg_factor ///< TODOCUMENT
		                                 ) {
			return {
				arg_coord.get_x() / arg_factor,
				arg_coord.get_y() / arg_factor,
				arg_coord.get_z() / arg_factor
			};
		}

		/// \brief TODOCUMENT
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::multiplicative<>
		///
		/// \relates coord
		inline constexpr coord operator/(const double &arg_factor, ///< TODOCUMENT
		                                 const coord  &arg_coord   ///< TODOCUMENT
		                                 ) {
			return {
				arg_coord.get_x() / arg_factor,
				arg_coord.get_y() / arg_factor,
				arg_coord.get_z() / arg_factor
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
		inline constexpr double squared_length(const coord &arg_coord ///< TODOCUMENT
		                                       ) {
			return (
				arg_coord.get_x() * arg_coord.get_x() +
				arg_coord.get_y() * arg_coord.get_y() +
				arg_coord.get_z() * arg_coord.get_z()
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline double length(const coord &arg_coord ///< TODOCUMENT
		                     ) {
			return sqrt( squared_length( arg_coord ) );
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline double distance_between_points(const coord &arg_coord1, ///< TODOCUMENT
		                                      const coord &arg_coord2  ///< TODOCUMENT
		                                      ) {
			// Equivalently, either :
			//   length(arg_coord2 - arg_coord1)
			// or
			//   sqrt(squared_distance_between_points(arg_coord1, arg_coord2));
			return length( arg_coord2 - arg_coord1 );
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline constexpr double squared_distance_between_points(const coord &arg_coord1, ///< TODOCUMENT
		                                                        const coord &arg_coord2  ///< TODOCUMENT
		                                                        ) {
			return squared_length( arg_coord2 - arg_coord1 );
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline constexpr bool operator==(const coord &arg_coord1, ///< TODOCUMENT
		                                 const coord &arg_coord2  ///< TODOCUMENT
		                                 ) {
			return (
				squared_distance_between_points( arg_coord1, arg_coord2 )
				<
				coord::TOLERANCE_FOR_COORD_CLOSENESS_CHECKS * coord::TOLERANCE_FOR_COORD_CLOSENESS_CHECKS
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline coord operator-(const coord &arg_coord ///< TODOCUMENT
		                       ) {
			return ( coord::ORIGIN_COORD - arg_coord );
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline bool not_zero(const coord &arg_coord ///< TODOCUMENT
		                     ) {
			if ( arg_coord.get_x() != 0.0 ) {
				return true;
			}
			if ( arg_coord.get_y() != 0.0 ) {
				return true;
			}
			if ( arg_coord.get_z() != 0.0 ) {
				return true;
			}
			return false;
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline double dot_product(const coord &arg_coord1, ///< TODOCUMENT
		                          const coord &arg_coord2  ///< TODOCUMENT
		                          ) {
			return (
				arg_coord1.get_x() * arg_coord2.get_x() +
				arg_coord1.get_y() * arg_coord2.get_y() +
				arg_coord1.get_z() * arg_coord2.get_z()
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline coord cross_product(const coord &arg_coord1, ///< TODOCUMENT
		                           const coord &arg_coord2  ///< TODOCUMENT
		                           ) {
			return coord(
				arg_coord1.get_y() * arg_coord2.get_z() - arg_coord1.get_z() * arg_coord2.get_y(),
				arg_coord1.get_z() * arg_coord2.get_x() - arg_coord1.get_x() * arg_coord2.get_z(),
				arg_coord1.get_x() * arg_coord2.get_y() - arg_coord1.get_y() * arg_coord2.get_x()
			);
		}

		/// \brief TODOCUMENT
		///
		/// \relates coord
		inline coord int_cast_copy(const coord &arg_coord ///< TODOCUMENT
		                           ) {
			return coord(
				boost::numeric_cast<int>( arg_coord.get_x() ),
				boost::numeric_cast<int>( arg_coord.get_y() ),
				boost::numeric_cast<int>( arg_coord.get_z() )
			);
			return arg_coord;
		}

		/// \brief Return the scalar component of the first specified coord in the direction of the second
		///
		/// \relates coord
		inline double scalar_component(const coord &arg_source_coord, ///< The coord of which to return the scalar component
		                               const coord &arg_dirn          ///< The direction in which the component should be taken (does not have to be normalised)
		                               ) {
			return dot_product( arg_source_coord, normalise_copy( arg_dirn ) );
		}

		/// \brief Return the parallel component of the first specified coord in the direction of the second
		///
		/// \relates coord
		inline coord parallel_component_copy(const coord &arg_source_coord, ///< The coord of which to return the parallel component
		                                     const coord &arg_dirn          ///< The direction in which the component should be taken (does not have to be normalised)
		                                     ) {
			return scalar_component( arg_source_coord, arg_dirn ) * normalise_copy( arg_dirn );
		}

		/// \brief Return the perpendicular component of the first specified coord in the direction perpendicular to the second
		///
		/// \relates coord
		inline coord perpendicular_component_copy(const coord &arg_source_coord, ///< The coord of which to return the perpendicular component
		                                          const coord &arg_dirn          ///< The direction perpendicular to which the component should be taken (does not have to be normalised)
		                                          ) {
			const coord  unit_dirn = normalise_copy( arg_dirn );
			const double factor    = dot_product   ( arg_source_coord, unit_dirn );
			return arg_source_coord - ( factor * unit_dirn );
		}


	} // namespace geom
	
	namespace common {
	
		/// \brief Specialisation of cath::common::read_from_ptree for coord
		template <>
		inline geom::coord read_from_ptree<geom::coord>(const boost::property_tree::ptree &arg_ptree ///< The ptree from which to read the coord
		                                                ) {
			return geom::coord_from_ptree( arg_ptree );
		}
	
	} // namespace common
} // namespace cath

#endif
