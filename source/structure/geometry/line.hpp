/// \file
/// \brief The eigen class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_LINE_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_GEOMETRY_LINE_H

#include "structure/geometry/coord.hpp"

namespace cath {
	namespace geom {

		/// \brief A simple line class
		///
		/// \todo This could be much improved, in particular
		///
		/// \todo Come relaxed constexpr support in all supported compilers,
		///       Improve this class (eg consider: setters, checking dirn has non-zero length)
		class line final {
		private:
			/// \brief A point on the line
			coord point_on_line;

			/// \brief A non-zero direction vector along the line (at present, this doesn't have to have length 1)
			coord dirn;

		public:
			line() = delete;

			/// \brief Ctor
			constexpr line(const coord &arg_point_on_line, ///< A point on the line
			               const coord &arg_dirn           ///< A non-zero direction vector along the line (at present, this doesn't have to have length 1)
			               ) : point_on_line{ arg_point_on_line },
			                   dirn         { arg_dirn          } {
			}

			/// \brief Getter for point_on_line
			const coord & get_point_on_line() const {
				return point_on_line;
			}

			/// \brief Getter for dirn
			const coord & get_dirn() const {
				return dirn;
			}
		};

		/// \brief Calculate the point on the specified line that is closest to the specified point
		inline coord closest_point_on_line_to_point(const line  &arg_line, ///< The line on which to find the point
		                                            const coord &arg_point ///< The point to get closest to
		                                            ) {
			return
				arg_line.get_point_on_line()
				+ parallel_component_copy(
					arg_point - arg_line.get_point_on_line(),
					arg_line.get_dirn()
				);
		}

	} // namespace chop
} // namespace cath

#endif
