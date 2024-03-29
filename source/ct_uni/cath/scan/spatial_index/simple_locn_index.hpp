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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SPATIAL_INDEX_SIMPLE_LOCN_INDEX_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SPATIAL_INDEX_SIMPLE_LOCN_INDEX_HPP

#include "cath/common/debug_numeric_cast.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"
#include "cath/structure/geometry/coord.hpp"

namespace cath::scan {

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
		constexpr simple_locn_index(const detail::view_base_type &prm_view_x, ///< TODOCUMENT
		                            const detail::view_base_type &prm_view_y, ///< TODOCUMENT
		                            const detail::view_base_type &prm_view_z, ///< TODOCUMENT
		                            const unsigned int           &prm_index   ///< TODOCUMENT
		                            ) : view_x { prm_view_x },
		                                view_y { prm_view_y },
		                                view_z { prm_view_z },
		                                index  { prm_index  } {
		}
	};

	/// \brief Convenience function to get the x component of the view in the specified simple_locn_index
	///
	/// \relates simple_locn_index
	inline constexpr const detail::view_base_type & get_view_x(const simple_locn_index &prm_locn_index ///< The simple_locn_index to query
	                                                           ) {
	        return prm_locn_index.view_x;
	}

	/// \brief Convenience function to get the y component of the view in the specified simple_locn_index
	///
	/// \relates simple_locn_index
	inline constexpr const detail::view_base_type & get_view_y(const simple_locn_index &prm_locn_index ///< The simple_locn_index to query
	                                                           ) {
	        return prm_locn_index.view_y;
	}

	/// \brief Convenience function to get the z component of the view in the specified simple_locn_index
	///
	/// \relates simple_locn_index
	inline constexpr const detail::view_base_type & get_view_z(const simple_locn_index &prm_locn_index ///< The simple_locn_index to query
	                                                           ) {
	        return prm_locn_index.view_z;
	}

	/// \brief Make a coord from the specified simple_locn_index
	///
	/// \relates simple_locn_index
	inline geom::coord get_coord(const simple_locn_index &prm_locn_index ///< The simple_locn_index to query
	                             ) {
		return {
			get_view_x( prm_locn_index ),
			get_view_y( prm_locn_index ),
			get_view_z( prm_locn_index )
		};
	}

	/// \brief Return whether the two specified full_hits are identical
	///
	/// \relates simple_locn_index
	inline bool operator==(const simple_locn_index &prm_locn_index_lhs, ///< The first  simple_locn_index to compare
	                       const simple_locn_index &prm_locn_index_rhs  ///< The second simple_locn_index to compare
	                       ) {
		return (
			get_view_x( prm_locn_index_lhs ) == get_view_x( prm_locn_index_rhs )
			&&
			get_view_y( prm_locn_index_lhs ) == get_view_y( prm_locn_index_rhs )
			&&
			get_view_z( prm_locn_index_lhs ) == get_view_z( prm_locn_index_rhs )
			&&
			prm_locn_index_lhs.index         == prm_locn_index_rhs.index
		);
	}

	/// \brief Get the squared distance between the two specified simple_locn_index values
	///
	/// \relates simple_locn_index
	inline detail::view_base_type get_squared_distance(const simple_locn_index &prm_locn_index_a, ///< The simple_locn_index to query
	                                                   const simple_locn_index &prm_locn_index_b  ///< The simple_locn_index to query
	                                                   ) {
		return (
			(
				( get_view_x( prm_locn_index_a ) - get_view_x( prm_locn_index_b ) )
				*
				( get_view_x( prm_locn_index_a ) - get_view_x( prm_locn_index_b ) )
			)
			+
			(
				( get_view_y( prm_locn_index_a ) - get_view_y( prm_locn_index_b ) )
				*
				( get_view_y( prm_locn_index_a ) - get_view_y( prm_locn_index_b ) )
			)
			+
			(
				( get_view_z( prm_locn_index_a ) - get_view_z( prm_locn_index_b ) )
				*
				( get_view_z( prm_locn_index_a ) - get_view_z( prm_locn_index_b ) )
			)
		);
	}

	/// \brief Return whether the two specified simple_locn_index values are within the specified distance
	///        (and associated squared distance)
	///
	/// \relates simple_locn_index
	inline bool are_within_distance_doub(const simple_locn_index &prm_locn_index_a,    ///< The simple_locn_index to query
	                                     const simple_locn_index &prm_locn_index_b,    ///< The simple_locn_index to query
	                                     const double            &prm_max_dist,        ///< The distance to which to compare
	                                     const double            &prm_max_squared_dist ///< The squared distance to which to compare
	                                     ) {
		const auto dist_x = debug_numeric_cast<double>( get_view_x( prm_locn_index_a ) ) - debug_numeric_cast<double>( get_view_x( prm_locn_index_b ) );
		if ( dist_x > prm_max_dist ) {
			return false;
		}
		const auto dist_y = debug_numeric_cast<double>( get_view_y( prm_locn_index_a ) ) - debug_numeric_cast<double>( get_view_y( prm_locn_index_b ) );
		if ( dist_y > prm_max_dist ) {
			return false;
		}
		const auto dist_z = debug_numeric_cast<double>( get_view_z( prm_locn_index_a ) ) - debug_numeric_cast<double>( get_view_z( prm_locn_index_b ) );
		if ( dist_z > prm_max_dist ) {
			return false;
		}
		return (
			( dist_x * dist_x )
			+
			( dist_y * dist_y )
			+
			( dist_z * dist_z )
		) < prm_max_squared_dist;
	}

	/// \brief TODOCUMENT
	///
	/// \relates simple_locn_index
	inline simple_locn_index make_simple_locn_index(const geom::coord  &prm_coord, ///< TODOCUMENT
	                                                const unsigned int &prm_index  ///< TODOCUMENT
	                                                ) {
		return {
			debug_numeric_cast< detail::view_base_type>( prm_coord.get_x() ),
			debug_numeric_cast< detail::view_base_type>( prm_coord.get_y() ),
			debug_numeric_cast< detail::view_base_type>( prm_coord.get_z() ),
			prm_index
		};
	}

	/// \brief TODOCUMENT
	struct simple_locn_crit final {

		/// \brief TODOCUMENT
		detail::view_base_type maximum_squared_distance;
	};

	inline std::string to_string(const simple_locn_index &prm_simple_locn_index
	                             ) {
		return "simple_locn_index[ ("
			+ ::std::to_string( prm_simple_locn_index.view_x )
			+ ", "
			+ ::std::to_string( prm_simple_locn_index.view_y )
			+ ", "
			+ ::std::to_string( prm_simple_locn_index.view_z )
			+ "), "
			+ ::std::to_string( prm_simple_locn_index.index )
			+ "]";
	}

	/// \brief TODOCUMENT
	inline detail::view_base_type get_maximum_distance(const simple_locn_crit &prm_crit ///< TODOCUMENT
	                                                   ) {
		return std::sqrt( prm_crit.maximum_squared_distance );
	}

} // namespace cath::scan

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SPATIAL_INDEX_SIMPLE_LOCN_INDEX_HPP
