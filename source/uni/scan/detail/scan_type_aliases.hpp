/// \file
/// \brief The scan type aliases header

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
///
/// The types defined in this header are key to scan speed so there is a greater
/// emphasis on compact types that can allow data to be packed more densely in vectors
///
/// \todo Split this file into scan's public type aliases and detail type aliases

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_SCAN_TYPE_ALIASES_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_SCAN_TYPE_ALIASES_HPP

#include <boost/geometry/arithmetic/arithmetic.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/optional/optional_fwd.hpp>
// #include <boost/units/base_units/information/byte.hpp>
#include <boost/units/systems/information.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "common/chrono/chrono_type_aliases.hpp"

namespace cath { namespace geom { template <typename T> class quat_rot_impl; } }
namespace cath { namespace geom { template <typename T> class angle; } }
namespace cath { namespace scan { namespace detail { class multi_struc_res_rep_pair; } } }
namespace cath { namespace scan { namespace detail { class multi_struc_res_rep_pair_list; } } }
namespace cath { namespace scan { namespace detail { class scan_structure_data; } } }
namespace cath { namespace scan { namespace detail { class single_struc_res_pair; } } }
namespace cath { namespace scan { namespace detail { class single_struc_res_pair_list; } } }
namespace cath { namespace scan { class quad_criteria; } }

namespace cath {
	namespace scan {
		/// \brief The type to be used for indexing structures / residues in res_pair code
		///
		/// It might be reasonable to consider uint16_t, which occupies two bytes. This would limit the maximum
		/// number of residues per structure to 65,536.
		///
		/// For comparison, mouse titin has 35,213 residues and here are some of the bigger pdbs
		/// and their rough numbers of residues:
		///     2cse 16,111
		///     1xi5 15,696
		///     3j8e 15,424
		///     1xi4 15,300
		///     3iyv 15,300
		///
		/// It's hard to imagine why anyone would want to compare whole structures of that magnitude.
		using index_type = uint32_t;

		/// \brief Type alias for a vector of the type to be used for indexing structures / residues in res_pair code
		using index_vec = std::vector<index_type>;

		namespace detail {
			/// \brief The type to be used for indexing the representative residues
			///
			/// Since there are typically substantially fewer representative residues
			/// than residues, this type can reasonably be implemented with fewer bytes
			using res_rep_index_type  = uint16_t;

			/// \brief Type alias for a pair of representative residue indices
			using rep_rep_pair        = std::pair<res_rep_index_type, res_rep_index_type>;

			/// \brief Type alias for an optional pair of representative residue indices
			using rep_rep_pair_opt    = boost::optional<rep_rep_pair>;

			/// \brief Type alias for a vector of scan_structure_data
			using scan_structure_data_vec      = std::vector<scan_structure_data>;

			/// \brief Type alias for a const_iterator for a vector of scan_structure_data
			using scan_structure_data_vec_citr = scan_structure_data_vec::const_iterator;

			/// \brief Type alias for a vector of multi_struc_res_rep_pair
			using multi_struc_res_rep_pair_vec      = std::vector<multi_struc_res_rep_pair>;

			/// \brief Type alias for a const_iterator for a vector of multi_struc_res_rep_pair
			using multi_struc_res_rep_pair_vec_citr = multi_struc_res_rep_pair_vec::const_iterator;

			/// \brief Type alias for a vector of multi_struc_res_rep_pair_list
			using multi_struc_res_rep_pair_list_vec = std::vector<multi_struc_res_rep_pair_list>;

			/// \brief Type alias for a vector of single_struc_res_pair
			using single_struc_res_pair_vec         = std::vector<single_struc_res_pair>;

			/// \brief Type alias for a const_iterator for a vector of single_struc_res_pair
			using single_struc_res_pair_vec_citr    = single_struc_res_pair_vec::const_iterator;

			/// \brief Type alias for a vector of single_struc_res_pair_list
			using single_struc_res_pair_list_vec = std::vector<single_struc_res_pair_list>;

			/// \brief Type alias template for pair of KEY and multi_struc_res_rep_pair_list
			template <typename KEY>
			using key_multi_res_pair_list_pair          = std::pair<KEY, multi_struc_res_rep_pair_list>;

			/// \brief Type alias template for vector of pair of KEY and multi_struc_res_rep_pair_list
			template <typename KEY>
			using key_multi_res_pair_list_pair_vec      = std::vector<key_multi_res_pair_list_pair<KEY>>;

			/// \brief Type alias template for the const_iterator of vector of pair of KEY and multi_struc_res_rep_pair_list
			template <typename KEY>
			using key_multi_res_pair_list_pair_vec_citr = typename key_multi_res_pair_list_pair_vec<KEY>::const_iterator;


			/// \brief The base type to be used in res_pair quaternion calculations about frames
			///        (meaning coordinate frames as defined by residues' atoms)
			using frame_quat_rot_type = float;

			/// \brief The quaternion type to be used in res_pair calculations about frames
			///        (meaning coordinate frames defined by residues' atoms)
			using frame_quat_rot      = geom::quat_rot_impl<frame_quat_rot_type>;

			/// \brief A vector of the type of quaternion to be used in res_pair calculations about frames
			///        (meaning coordinate frames defined by residues' atoms)
			using frame_quat_rot_vec  = std::vector<frame_quat_rot>;

			/// \brief The base type to be used in angle<> for representing angles in res_pair code
			using angle_base_type     = float;
//			using angle_base_type     = double;

			/// \brief The full angle<> type to be used for representing angles in res_pair code
			using angle_type          = geom::angle<angle_base_type>;

			/// \brief A vector of the angle_type
			using angle_type_vec      = std::vector<angle_type>;

			/// \brief The base type to be used in coords for representing views in res_pair
			using view_base_type      = float;
//			using view_base_type      = double;

			/// \brief The type used to represent the view
//			using view_type           = geom::coord;
			using view_type           = boost::geometry::model::point<view_base_type, 3, boost::geometry::cs::cartesian>;

			inline view_type operator-(view_type        prm_view_a,
			                           const view_type &prm_view_b
			                           ) {
				boost::geometry::subtract_point( prm_view_a, prm_view_b );
				return prm_view_a;
			}

			/// \brief The type for a vector of view_types
			using view_type_vec       = std::vector<view_type>;

			/// \brief The default numeric type used to multiply types
			using multiplier_type       = float;

			/// \brief The default signed int type used for indices for x,y or z dimensions of the view
			using key_view_index_type   = int16_t;

			/// \brief The default unsigned int type used for indices for the orientation
			using key_orient_index_type = uint8_t;

			/// \brief The default unsigned int type used for indices for from/to phi/psi angles
			using key_angle_index_type  = uint8_t;
		} // namespace detail

		/// \brief TODOCUMENT
		using info_value        = size_t;

		/// \brief TODOCUMENT
		using info_quantity     = boost::units::quantity<boost::units::information::info, info_value>;

		/// \brief TODOCUMENT
		using durn_mem_pair     = std::pair<hrc_duration, info_quantity>;

		// /// \brief TODOCUMENT
		// using durn_mem_pair_opt = boost::optional<durn_mem_pair>;

		/// \brief TODOCUMENT
		using hrc_duration_opt  = boost::optional<hrc_duration>;

		/// \brief Type alias for a vector of quad_criteria
		using quad_criteria_vec = std::vector<quad_criteria>;
	} // namespace scan
} // namespace cath

#endif
