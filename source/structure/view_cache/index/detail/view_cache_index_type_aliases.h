/// \file
/// \brief The view_cache_index class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef VIEW_CACHE_INDEX_TYPE_ALIASES_H_INCLUDED
#define VIEW_CACHE_INDEX_TYPE_ALIASES_H_INCLUDED

#include <boost/geometry/geometries/geometries.hpp> // ***** TEMPORARY? *****
#include <boost/geometry/geometries/point.hpp> // ***** TEMPORARY? *****

namespace cath { namespace geom { template <typename T> class quat_rot_impl; } }
namespace cath { namespace geom { class coord; } }

namespace cath {
	namespace index {
		class view_cache_index_entry;

		/// \brief Convenience type alias for a vector of view_cache_index_entry
		using view_cache_index_entry_vec = std::vector<view_cache_index_entry>;

		/// \brief Convenience type alias for a vector of view_cache_index_entry
		using vcie_vcie_vec_pair         = std::pair<view_cache_index_entry_vec, view_cache_index_entry_vec>;

		/// \brief Convenience type alias for a vector of view_cache_index_entry
		using vcie_vcie_pair             = std::pair<view_cache_index_entry, view_cache_index_entry>;

		/// \brief Convenience type alias for a vector of view_cache_index_entry
		using vcie_vcie_pair_vec         = std::vector<vcie_vcie_pair>;

		namespace detail {

			/// \brief The base type to be used in view_cache_index quaternion calculations about frames
			///        (meaning coordinate frames as defined by residues' atoms)
			using frame_quat_rot_type = float;

			/// \brief The quaternion type to be used in view_cache_index calculations about frames
			///        (meaning coordinate frames defined by residues' atoms)
			using frame_quat_rot      = geom::quat_rot_impl<frame_quat_rot_type>;

			/// \brief A vector of the type of quaternion to be used in view_cache_index calculations about frames
			///        (meaning coordinate frames defined by residues' atoms)
			using frame_quat_rot_vec  = std::vector<frame_quat_rot>;

			/// \brief The type to be used for representing indices in view_cache_index code
			///
			/// unsigned int (which typically occupies 32 bits/4 bytes on Linux) is used here
			/// in preference to size_t (which typically occupies 64 bits/8 bytes on 64-bit Linux)
			/// because speed (and thus memory compactness) is very important in view_cache_index code.
			///
			/// \todo Consider uint16_t, which occupies two bytes but which would limit the maximum
			///       number of residues per structure to 65,536.
			///       For comparison, mouse titin has 35,213 residues
			///
			/// Here are some of the bigger pdbs and their numbers of residues:
			///     2cse 16111
			///     1xi5 15696
			///     3j8e 15424
			///     1xi4 15300
			///     3iyv 15300
			using index_type          = unsigned int;

			/// \brief The base type to be used in angle<> for representing angles in view_cache_index code
			using angle_base_type     = float;
//			using angle_base_type     = double;

			/// \brief The full angle<> type to be used for representing angles in view_cache_index code
			using angle_type          = geom::angle<angle_base_type>;

			/// \brief The base type to be used in coords for representing views in view_cache_index
			using view_base_type      = float;
//			using view_base_type      = double;

			/// \brief The type used to represent the view
//			using view_type           = geom::coord;
			using view_type           = boost::geometry::model::point<view_base_type, 3, boost::geometry::cs::cartesian>;

			/// \brief The default numeric type used to multiply types
			using multiplier_type     = float;
		}
	}
}

#endif
