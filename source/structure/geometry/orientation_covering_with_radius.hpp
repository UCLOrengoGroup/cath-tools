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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_ORIENTATION_COVERING_WITH_RADIUS_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_ORIENTATION_COVERING_WITH_RADIUS_H

// #include <boost/range/adaptor/transformed.hpp>
// #include <boost/range/adaptor/filtered.hpp>
// #include <boost/range/adaptor/transformed.hpp>
// #include <boost/range/algorithm/min_element.hpp>
// #include <boost/algorithm/string/classification.hpp>
// #include <boost/algorithm/string/predicate.hpp>
// #include <boost/algorithm/string/trim.hpp>
// #include <boost/lexical_cast.hpp>
// #include <boost/range/irange.hpp>

// #include "common/algorithm/copy_build.hpp"
// #include "common/algorithm/sort_uniq_copy.hpp"
// #include "common/algorithm/transform_build.hpp"
// #include "common/boost_addenda/range/range_concept_type_aliases.hpp"
// #include "common/boost_addenda/string_algorithm/split_build.hpp"
// #include "common/debug_numeric_cast.hpp"
// #include "common/file/open_fstream.hpp"
// #include "common/size_t_literal.hpp"
// #include "structure/geometry/quat_rot.hpp"
// #include "exception/invalid_argument_exception.hpp"
// #include "exception/runtime_error_exception.hpp"
#include "structure/geometry/orientation_covering.hpp"
// #include "structure/structure_type_aliases.hpp"

// // #include <algorithm>
// // #include <cmath>
// #include <iostream> // ***** TEMPORARY *****
// #include <fstream>

// using namespace cath::common::literals;

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
		inline orientation_covering_with_radius<T>::orientation_covering_with_radius(const angle<T> &arg_search_radius ///< TODOCUMENT
		                                                                             ) : search_radius   ( arg_search_radius ),
		                                                                                 cell_neighbours ( calc_neighbours( the_covering, search_radius ) ) {
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline uint8_t orientation_covering_with_radius<T>::get_closest_neighbour(const quat_rot_impl<T> &arg_orientation ///< TODOCUMENT
		                                                                          ) const {
			return get_closest_neighbour( the_covering, arg_orientation );
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline uint8_vec orientation_covering_with_radius<T>::get_closest_neighbours(const quat_rot_impl<T> &arg_orientation ///< TODOCUMENT
		                                                                             ) const {
			return get_closest_neighbours( the_covering, arg_orientation, cell_neighbours, search_radius );
		}

		/// \brief TODOCUMENT
		template <typename T>
		inline const angle<T> & orientation_covering_with_radius<T>::get_search_radius() const {
			return search_radius;
		}


	} // namespace geom
} // namespace cath

#endif
