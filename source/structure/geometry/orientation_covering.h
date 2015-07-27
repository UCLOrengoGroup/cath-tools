/// \file
/// \brief The orientation_covering class header

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

#ifndef ORIENTATION_COVERING_H_INCLUDED
#define ORIENTATION_COVERING_H_INCLUDED

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/copy_build.h"
#include "common/algorithm/sort_uniq_copy.h"
#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/range/range_concept_type_aliases.h"
#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/debug_numeric_cast.h"
#include "common/file/open_fstream.h"
#include "common/size_t_literal.h"
#include "structure/geometry/quat_rot.h"
#include "exception/invalid_argument_exception.h"
#include "exception/runtime_error_exception.h"
#include "structure/structure_type_aliases.h"

// #include <algorithm>
// #include <cmath>
#include <iostream> // ***** TEMPORARY *****
#include <fstream>


using namespace cath::common::literals;

namespace cath {
	namespace geom {

		/// \brief A set of orientations that cover orientation space
		///
		/// The data sets aim to provide a nearly optimal covering of orientation
		/// space (ie using as few orientations as possible such that all of orientation
		/// space is within some specific angle of one of those orientations)
		///
		/// (The data's from http://charles.karney.info/orientation/ - see that URL for more info)
		template <typename T>
		class orientation_covering_impl final {
		private:
			/// \brief The set of orientations that form the covering set
			quat_rot_vec<T> orientations;

			/// \brief TODOCUMENT
			static boost::filesystem::path get_orientation_filename();

			/// \brief Load a list of orientations from a file
			///
			/// This is currently used to initialise the orientations on construction
			static quat_rot_vec<T> orientations_from_file(const boost::filesystem::path &);

		public:
			/// \brief TODOCUMENT
			using const_iterator = common::range_const_iterator_t<quat_rot_vec<T>>;

			orientation_covering_impl();

			bool empty() const;
			size_t size() const;
			const quat_rot_impl<T> & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		/// \brief TODOCUMENT
		template <typename T>
		quat_rot_vec<T> get_c600v() {
			// grep -vP '^#' c600v.quat | awk '$5 > 0 {printf( "{ %12s, %12s, %12s, %12s  },\n", $1, $2, $3, $4 ); }'
			return {
				quat_rot_impl<T>{ static_cast<T>(  1.000000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  1.000000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  1.000000000 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  1.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  0.500000000 ), static_cast<T>( -0.309016994 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.000000000 ), static_cast<T>( -0.500000000 ), static_cast<T>( -0.309016994 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>( -0.309016994 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.000000000 ), static_cast<T>( -0.309016994 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.000000000 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.500000000 ), static_cast<T>( -0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.000000000 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.500000000 ), static_cast<T>( -0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.500000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.500000000 ), static_cast<T>( -0.000000000 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.309016994 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.309016994 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.309016994 ), static_cast<T>(  0.000000000 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.309016994 ), static_cast<T>( -0.000000000 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>( -0.000000000 ), static_cast<T>( -0.309016994 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>(  0.000000000 ), static_cast<T>( -0.309016994 ), static_cast<T>(  0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.309016994 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.309016994 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.309016994 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.309016994 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.000000000 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.309016994 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.000000000 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>( -0.500000000 ), static_cast<T>( -0.000000000 ), static_cast<T>( -0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.000000000 ), static_cast<T>( -0.309016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.500000000 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>( -0.309016994 ), static_cast<T>(  0.000000000 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.500000000 ), static_cast<T>( -0.309016994 ), static_cast<T>( -0.000000000 ), static_cast<T>(  0.809016994 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>( -0.309016994 ), static_cast<T>( -0.500000000 ), static_cast<T>( -0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>( -0.309016994 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.809016994 ), static_cast<T>(  0.309016994 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.500000000 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.309016994 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.500000000 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.309016994 ), static_cast<T>( -0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.309016994 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.309016994 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.309016994 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.309016994 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.809016994 ), static_cast<T>(  0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>( -0.309016994 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.809016994 ), static_cast<T>( -0.000000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>( -0.500000000 ), static_cast<T>( -0.500000000 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>(  0.500000000 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>(  0.500000000 ), static_cast<T>(  0.500000000 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>(  0.500000000 ), static_cast<T>( -0.500000000 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.500000000 ), static_cast<T>( -0.500000000 ) },
				quat_rot_impl<T>{ static_cast<T>(  0.500000000 ), static_cast<T>( -0.500000000 ), static_cast<T>( -0.500000000 ), static_cast<T>(  0.500000000 ) }
			};
		}

// 		/// \brief TODOCUMENT
// 		template <typename T>
// 		boost::filesystem::path orientation_covering_impl<T>::get_orientation_filename() {
// //			return { "/tmp/c48u1.quat"   }; // num:  24; max_angle: 62.80
// 			return { "/tmp/c600v.quat"   }; // num:  60; max_angle: 44.48
// //			return { "/tmp/c48n9.quat"   }; // num: 216; max_angle: 36.47
// //			return { "/tmp/c48u9.quat"   }; // num: 216; max_angle: 38.45
// //			return { "/tmp/c600vc.quat"  }; // num: 360; max_angle: 27.78
// //			return { "/tmp/c48u27.quat"  }; // num: 648; max_angle: 20.83
// //			return { "/tmp/c600vec.quat" }; // num: 720; max_angle: 22.25
// 		}

		/// \brief TODOCUMENT
		template <typename T>
		quat_rot_vec<T> orientation_covering_impl<T>::orientations_from_file(const boost::filesystem::path &arg_filename ///< TODOCUMENT
		                                                                     ) {
			// Open an ifstream for the file
			std::ifstream orentation_ifstream;
			std::cerr << "About to attempt to open file " << arg_filename.string() << std::endl;
			common::open_ifstream( orentation_ifstream, arg_filename );

			/// Populate the orientations, one per line
			quat_rot_vec<T> orientations;
			std::string line_string;
			while ( std::getline( orentation_ifstream, line_string ) ) {
				// Trim off any leading/trailing whitespace
				boost::algorithm::trim( line_string );

				// If this is a comment line, then skip to the next line
				if ( boost::algorithm::starts_with( line_string, "#" ) ) {
					continue;
				}
				// If this is a format line then check the format and then skip this line and the next
				if ( boost::algorithm::starts_with( line_string, "format " ) ) {
					if ( ! boost::algorithm::starts_with( line_string, "format quaternion" ) ) {

					}
					std::getline( orentation_ifstream, line_string );
					continue;
				}

				// Otherwise this should be a data line so split it into parts on whitespace
				const str_vec quat_str_parts = common::split_build<str_vec>(
					line_string,
					boost::algorithm::is_space(),
					boost::algorithm::token_compress_on
				);

				// If there aren't five parts to the line then throw an error
				if ( quat_str_parts.size() != 5 ) {
					constexpr size_t ERROR_LINE_LENGTH = 1000;
					BOOST_THROW_EXCEPTION(common::runtime_error_exception(
						"Whilst attempting to read orientation covering file "
						+ arg_filename.string()
						+ ", unable to parse five parts from line: \""
						+ line_string.substr( 0, ERROR_LINE_LENGTH - 1 ) + ( line_string.length() > ERROR_LINE_LENGTH ? "[...]" : "" )
						+ "\""
					));
				}

				// Build a quaternion from the first four numbers
				const auto quat_parts = common::transform_build<std::vector<T>>( quat_str_parts, &boost::lexical_cast<T, std::string> );
				orientations.emplace_back(
					quat_parts[ 0 ],
					quat_parts[ 1 ],
					quat_parts[ 2 ],
					quat_parts[ 3 ]
				);

//				std::cerr << "Read line "   << line_string
//				          << ", real is "   << real  ( orientations.back() )
//				          << ", unreal is " << unreal( orientations.back() )
//				          << ", sup is "    << sup   ( orientations.back() )
//				          << ", l1 is "     << l1    ( orientations.back() )
//				          << ", abs is "    << abs   ( orientations.back() )
//				          << ", norm is "   << norm  ( orientations.back() )
//				          << "\n";
			}

			// Close the ifstream to check that there are no errors
			orentation_ifstream.close();

			// Return the loaded orientations
			return orientations;
		}

		/// \brief TODOCUMENT
		template <typename T>
		orientation_covering_impl<T>::orientation_covering_impl() : orientations( get_c600v<T>() ) {
		}

		/// \brief TODOCUMENT
		template <typename T>
		bool orientation_covering_impl<T>::empty() const {
			return orientations.empty();
		}

		/// \brief TODOCUMENT
		template <typename T>
		size_t orientation_covering_impl<T>::size() const {
			return orientations.size();
		}

		/// \brief TODOCUMENT
		template <typename T>
		const quat_rot_impl<T> & orientation_covering_impl<T>::operator[](const size_t &arg_index ///< TODOCUMENT
		                                                                  ) const {
			return orientations[ arg_index ];
		}

		/// \brief TODOCUMENT
		template <typename T>
		auto orientation_covering_impl<T>::begin() const -> const_iterator {
			return common::cbegin( orientations );
		}

		/// \brief TODOCUMENT
		template <typename T>
		auto orientation_covering_impl<T>::end() const -> const_iterator {
			return common::cend( orientations );
		}

		/// \brief TODOCUMENT
		///
		/// How to determine whether two adjacent covering orientations should be considered
		/// neighbours (for a given radius)?
		///
		/// This is made trickier (than for, eg, the 3d grid for the view) because:
		///  * different pairs of cells may have different distances (ie angles) between their cell-centres
		///  * a pair of cells that do contain points within arg_radius of each other might have a third
		///    cell directly in between their mid-points.
		///
		///
		/// Currently misses one every ~7142.85714286 of random pairs about 21.5 - 23.5 degrees apart
		/// under a radius slightly larger than their distance apart.
		///
		/// and because the line between two cell-centres may cross another cell but only for some angle
		/// less than the radius (meaning that there may be two points in the two cells that are within radius of
		/// each other.
		///
		///
		/// NOTE: First draft of this attempted to find whether each of the two centre points
		///       were within 2*radius of each other. This isn't adequate because two cells
		///       should be considered neighbours under the weaker condition that there exist
		///       two points with radius of each other such that each is closer to one of the
		///       two query cell-centres than to any other.
		///
		/// Note: this makes a few mathematical assumptions that I haven't completely checked
		///       so it's important to test that this calculates correct results...
		///
		/// \todo Test this by scanning the all-against-all anchor pairs and then for each, check:
		///        * the half-way quaternion is the same angle to both of the pair
		///        * the pair's in neighbours iff the half-way quaternion is within arg_search_radius of both
		///
		/// \todo Once the above tests are in place, try changing the code to calculate the
		///       distance_1_of_angle(2.0 * arg_search_radius) at the start and then compare
		///       the distance_1_of_quat_rot( query_orientn, match_orientn ) to that value
		template <typename T, typename A>
		size_vec_vec calc_neighbours(const orientation_covering_impl<T> &arg_orientations, ///< TODOCUMENT
		                             const geom::angle<A>               &arg_search_radius ///< TODOCUMENT
		                             ) {
			const auto orientation_index_range = boost::irange( 0_z, arg_orientations.size() );

			return common::transform_build<size_vec_vec>(
				orientation_index_range,
				[&] (const size_t &query_orientn_idx) {
					return common::copy_build<size_vec>(
						orientation_index_range
						| boost::adaptors::filtered(
							[&] (const size_t &match_orientn_idx) {
								const auto &query_orientn   = arg_orientations[ query_orientn_idx ];
								const auto &match_orientn   = arg_orientations[ match_orientn_idx ];
								const auto  angle_between   = angle_between_quat_rots( query_orientn, match_orientn );
								if ( angle_between <= arg_search_radius ) {
									return true;
								}
								const auto  half_between    = angle_between     / 2.0;
								const auto  half_radius     = arg_search_radius / 2.0;
								const auto  first_angle     = half_between - half_radius;
								const auto  second_angle    = half_between + half_radius;;
								const auto  first_point     = from_first_toward_second_at_angle( query_orientn, match_orientn, first_angle  );
								const auto  second_point    = from_first_toward_second_at_angle( query_orientn, match_orientn, second_angle );
								const bool  first_point_in  = get_closest_neighbour( arg_orientations, first_point  ) == query_orientn_idx;
								const bool  second_point_in = get_closest_neighbour( arg_orientations, second_point ) == match_orientn_idx;
								return ( first_point_in && second_point_in );
							}
						)
					);
				}
			);
		}

		/// \brief TODOCUMENT
		template <typename T, typename A>
		auto get_closest_neighbours(const orientation_covering_impl<T> &arg_orientations, ///< TODOCUMENT
		                            const quat_rot_impl<T>             &arg_orientation,  ///< TODOCUMENT
		                            const size_vec_vec                 &arg_neighbours,   ///< TODOCUMENT
		                            const geom::angle<A>               &arg_search_radius ///< TODOCUMENT
		                            ) {
			const size_t  closest_neighbour_idx = get_closest_neighbour( arg_orientations, arg_orientation );
			const auto   &closest_neighbour     = arg_orientations[ closest_neighbour_idx ];
//			std::cerr << "The original number of neighbours is " << arg_neighbours[ closest_neighbour_idx ].size() << "\n";
			return common::copy_build<size_vec>(
				arg_neighbours[ closest_neighbour_idx ]
					| boost::adaptors::filtered(
						[&] (const size_t &x) {
							if ( x == closest_neighbour_idx ) {
								return true;
							}
							const auto &neighbour = arg_orientations[ x ];
//							std::cerr << "considering " << neighbour;
							if ( angle_between_quat_rots( arg_orientation, neighbour ) <= arg_search_radius ) {
//								std::cerr << " - accepted immediately\n";
								return true;
							}
							const auto point1 = from_first_toward_second_at_angle( arg_orientation, neighbour, arg_search_radius );
							const bool result = ( distance_1_between_quat_rots( point1, neighbour ) < distance_1_between_quat_rots( point1, closest_neighbour ) );
//							std::cerr << " - considering " << point << ", which is " << angle_between_quat_rots( point, arg_orientation ) << " from target - result" << std::boolalpha << result << "\n";
							if ( result ) {
								return true;
							}
//							return ( distance_1_between_quat_rots( point, neighbour ) < distance_1_between_quat_rots( point, closest_neighbour ) );

							const auto last_ditch_angle = std::min( arg_search_radius, angle_between_quat_rots(  closest_neighbour, neighbour ) );
							const auto A_to_B_by_r      = from_first_toward_second_at_angle(  closest_neighbour, neighbour, last_ditch_angle );
							const auto dirn_A_to_b_by_r = rotation_between_rotations( closest_neighbour, A_to_B_by_r );
							const auto point2           = quat_rot_impl<T>{ dirn_A_to_b_by_r * arg_orientation };
							return ( distance_1_between_quat_rots( point2, neighbour ) < distance_1_between_quat_rots( point2, closest_neighbour ) );
						}
					)
			);
		}

		/// \brief TODOCUMENT
		///
		/// This could likely be made more efficient, particularly if a particular quaterion
		/// covering is settled on.
		template <typename T>
		size_t get_closest_neighbour(const orientation_covering_impl<T> &arg_orientations, ///< TODOCUMENT
		                             const quat_rot_impl<T>             &arg_orientation   ///< TODOCUMENT
		                             ) {
			// Check that the orientations aren't empty
			if ( arg_orientations.empty() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to get closest neighbour from empty orientation covering"));
			}

			// Search for the orientation_covering entry closest to arg_orientation
			const auto min_itr = boost::range::min_element(
				arg_orientations,
				[&] (const quat_rot_impl<T> &x, const quat_rot_impl<T> &y) {
					return ( distance_1_between_quat_rots( x, arg_orientation ) < distance_1_between_quat_rots( y, arg_orientation ) );
				}
			);
//			std::cerr << "Found " << *min_itr << "\n";
			
			// Check that some entry was chosen
			assert( min_itr != common::cend( arg_orientations ) );

			// Return the index of the element found
			return debug_numeric_cast<size_t>( distance( common::cbegin( arg_orientations ), min_itr ) );
		}

		


	}
}

#endif
