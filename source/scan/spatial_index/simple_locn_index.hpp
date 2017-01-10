/// \file
/// \brief The simple_locn_index class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_SPATIAL_INDEX_SIMPLE_LOCN_INDEX_H
#define _CATH_TOOLS_SOURCE_SCAN_SPATIAL_INDEX_SIMPLE_LOCN_INDEX_H

#include "common/debug_numeric_cast.hpp"
#include "scan/detail/scan_type_aliases.hpp"
#include "structure/geometry/coord.hpp"

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		struct simple_locn_index final {
			using view_t = detail::view_base_type;

			/// \brief TODOCUMENT
			detail::view_base_type view_x;

			/// \brief TODOCUMENT
			detail::view_base_type view_y;

			/// \brief TODOCUMENT
			detail::view_base_type view_z;

			/// \brief TODOCUMENT
			unsigned int index;

			/// \brief TODOCUMENT
			constexpr simple_locn_index(const detail::view_base_type &arg_view_x, ///< TODOCUMENT
			                            const detail::view_base_type &arg_view_y, ///< TODOCUMENT
			                            const detail::view_base_type &arg_view_z, ///< TODOCUMENT
			                            const unsigned int           &arg_index   ///< TODOCUMENT
			                            ) : view_x { arg_view_x },
			                                view_y { arg_view_y },
			                                view_z { arg_view_z },
			                                index  { arg_index  } {
			}
		};

		/// \brief Convenience function to get the x component of the view in the specified simple_locn_index
		///
		/// \relates simple_locn_index
		inline constexpr const detail::view_base_type & get_view_x(const simple_locn_index &arg_locn_index ///< The simple_locn_index to query
		                                                           ) {
		        return arg_locn_index.view_x;
		}

		/// \brief Convenience function to get the y component of the view in the specified simple_locn_index
		///
		/// \relates simple_locn_index
		inline constexpr const detail::view_base_type & get_view_y(const simple_locn_index &arg_locn_index ///< The simple_locn_index to query
		                                                           ) {
		        return arg_locn_index.view_y;
		}

		/// \brief Convenience function to get the z component of the view in the specified simple_locn_index
		///
		/// \relates simple_locn_index
		inline constexpr const detail::view_base_type & get_view_z(const simple_locn_index &arg_locn_index ///< The simple_locn_index to query
		                                                           ) {
		        return arg_locn_index.view_z;
		}

		inline detail::view_base_type get_squared_distance(const simple_locn_index &arg_locn_index_a, ///< The simple_locn_index to query
		                                                   const simple_locn_index &arg_locn_index_b  ///< The simple_locn_index to query
		                                                   ) {
			return (
				(
					( get_view_x( arg_locn_index_a ) - get_view_x( arg_locn_index_b ) )
					*
					( get_view_x( arg_locn_index_a ) - get_view_x( arg_locn_index_b ) )
				)
				+
				(
					( get_view_y( arg_locn_index_a ) - get_view_y( arg_locn_index_b ) )
					*
					( get_view_y( arg_locn_index_a ) - get_view_y( arg_locn_index_b ) )
				)
				+
				(
					( get_view_z( arg_locn_index_a ) - get_view_z( arg_locn_index_b ) )
					*
					( get_view_z( arg_locn_index_a ) - get_view_z( arg_locn_index_b ) )
				)
			);
		}

		inline bool are_within_distance(const simple_locn_index &arg_locn_index_a, ///< The simple_locn_index to query
		                                const simple_locn_index &arg_locn_index_b, ///< The simple_locn_index to query
		                                const float &arg_max_dist,                 ///< TODOCUMENT
		                                const float &arg_max_squared_dist          ///< TODOCUMENT
		                                ) {
			const auto dist_x = get_view_x( arg_locn_index_a ) - get_view_x( arg_locn_index_b );
			if ( dist_x > arg_max_dist ) {
				return false;
			}
			const auto dist_y = get_view_y( arg_locn_index_a ) - get_view_y( arg_locn_index_b );
			if ( dist_y > arg_max_dist ) {
				return false;
			}
			const auto dist_z = get_view_z( arg_locn_index_a ) - get_view_z( arg_locn_index_b );
			if ( dist_z > arg_max_dist ) {
				return false;
			}
			return (
				( dist_x * dist_x )
				+
				( dist_y * dist_y )
				+
				( dist_z * dist_z )
			) < arg_max_squared_dist;
		}

		/// \brief TODOCUMENT
		inline simple_locn_index make_simple_locn_index(const geom::coord  &arg_coord, ///< TODOCUMENT
		                                                const unsigned int &arg_index  ///< TODOCUMENT
		                                                ) {
			return {
				debug_numeric_cast< detail::view_base_type>( arg_coord.get_x() ),
				debug_numeric_cast< detail::view_base_type>( arg_coord.get_y() ),
				debug_numeric_cast< detail::view_base_type>( arg_coord.get_z() ),
				arg_index
			};
		}

		/// \brief TODOCUMENT
		struct simple_locn_crit final {

			/// \brief TODOCUMENT
			detail::view_base_type maximum_squared_distance;
		};

		inline std::string to_string(const simple_locn_index &arg_simple_locn_index
		                             ) {
			return "simple_locn_index[ ("
				+ ::std::to_string( arg_simple_locn_index.view_x )
				+ ", "
				+ ::std::to_string( arg_simple_locn_index.view_y )
				+ ", "
				+ ::std::to_string( arg_simple_locn_index.view_z )
				+ "), "
				+ ::std::to_string( arg_simple_locn_index.index )
				+ "]";
		}

		/// \brief TODOCUMENT
		inline detail::view_base_type get_maximum_distance(const simple_locn_crit &arg_crit ///< TODOCUMENT
		                                                   ) {
			return std::sqrt( arg_crit.maximum_squared_distance );
		}

	} // namespace scan
} // namespace cath

#endif
