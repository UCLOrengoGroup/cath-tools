/// \file
/// \brief The orientation_covering with radius class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ORIENTATION_COVERING_WITH_RADIUS_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ORIENTATION_COVERING_WITH_RADIUS_HPP

// #include <boost/range/adaptor/transformed.hpp>
// #include <boost/range/adaptor/filtered.hpp>
// #include <boost/range/adaptor/transformed.hpp>
// #include <boost/range/algorithm/min_element.hpp>
// #include <boost/algorithm/string/classification.hpp>
// #include <boost/algorithm/string/predicate.hpp>
// #include <boost/algorithm/string/trim.hpp>
// #include <boost/lexical_cast.hpp>

#include "cath/structure/geometry/orientation_covering.hpp"
// #include "cath/common/algorithm/copy_build.hpp"
// #include "cath/common/algorithm/sort_uniq_copy.hpp"
// #include "cath/common/algorithm/transform_build.hpp"
// #include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
// #include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
// #include "cath/common/debug_numeric_cast.hpp"
// #include "cath/common/exception/invalid_argument_exception.hpp"
// #include "cath/common/exception/runtime_error_exception.hpp"
// #include "cath/common/file/open_fstream.hpp"
// #include "cath/structure/geometry/quat_rot.hpp"
// #include "cath/structure/structure_type_aliases.hpp"

// // #include <algorithm>
// // #include <cmath>
// #include <iostream> // ***** TEMPORARY *****
// #include <fstream>


namespace cath {
	namespace geom {

		/// \brief TODOCUMENT
		using uint8_vec     = std::vector<uint8_t>;

		/// \brief TODOCUMENT
		using uint8_vec_vec = std::vector<uint8_vec>;

		/// \brief TODOCUMENT
		template <typename T>
		class orientation_covering_with_radius final {
		private:
			/// \brief The covering set
			orientation_covering_impl<T> the_covering;

			/// \brief TODOCUMENT
			angle<T> search_radius;

			/// \brief TODOCUMENT
			uint8_vec_vec cell_neighbours;

		public:
			explicit orientation_covering_with_radius(const angle<T> &);

			uint8_t get_closest_neighbour(const quat_rot_impl<T> &) const;
			uint8_vec get_closest_neighbours(const quat_rot_impl<T> &) const;

			const angle<T> & get_search_radius() const;
		};

		/// \brief TODOCUMENT
		template <typename T>
		inline orientation_covering_with_radius<T>::orientation_covering_with_radius(const angle<T> &prm_search_radius ///< TODOCUMENT
		                                                                             ) : search_radius   ( prm_search_radius ),
		                                                                                 cell_neighbours ( calc_neighbours( the_covering, search_radius ) ) {
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline uint8_t orientation_covering_with_radius<T>::get_closest_neighbour(const quat_rot_impl<T> &prm_orientation ///< TODOCUMENT
		                                                                          ) const {
			return get_closest_neighbour( the_covering, prm_orientation );
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline uint8_vec orientation_covering_with_radius<T>::get_closest_neighbours(const quat_rot_impl<T> &prm_orientation ///< TODOCUMENT
		                                                                             ) const {
			return get_closest_neighbours( the_covering, prm_orientation, cell_neighbours, search_radius );
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline const angle<T> & orientation_covering_with_radius<T>::get_search_radius() const {
			return search_radius;
		}


	} // namespace geom
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_GEOMETRY_ORIENTATION_COVERING_WITH_RADIUS_HPP
