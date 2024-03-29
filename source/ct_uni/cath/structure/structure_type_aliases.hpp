/// \file
/// \brief The structure type_aliases header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_STRUCTURE_TYPE_ALIASES_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_STRUCTURE_TYPE_ALIASES_HPP

#include "cath/common/type_aliases.hpp"
#include "cath/structure/geometry/coord_linkage.hpp"

#include <set>
#include <vector>

// clang-format off
namespace cath { class amino_acid; }
namespace cath { class chain_label; }
namespace cath { class protein; }
namespace cath { class residue; }
namespace cath { class residue_id; }
namespace cath { class residue_name; }
namespace cath { class sec_struc; }
namespace cath { class sec_struc_planar_angles; }
namespace cath { enum class sec_struc_type : char; }
namespace cath::geom { class coord; }
namespace cath::geom { class coord_list; }
namespace cath::geom { class rotation; }
namespace cath::geom { template <typename T> class angle; }
namespace cath::geom { template <typename T> class quat_rot_impl; }
namespace cath::index { class view_cache; }
// clang-format on

namespace cath {
	using sec_struc_type_vec = std::vector<sec_struc_type>;
} // namespace cath

namespace cath {
	namespace index {

		namespace detail {
			class vcie_match_criteria;

			/// \brief TODOCUMENT
			using vcie_match_criteria_vec = std::vector<vcie_match_criteria>;
		} // namespace detail

		namespace filter {
			class filter_vs_full_score;

			/// \brief Type alias for vector of filter_vs_full_scores
			using filter_vs_full_score_vec      = std::vector<filter_vs_full_score>;

			/// \brief Type alias for filter_vs_full_score_vec's iterator
			using filter_vs_full_score_vec_itr  = filter_vs_full_score_vec::iterator;

			/// \brief Type alias for filter_vs_full_score_vec's const_iterator
			using filter_vs_full_score_vec_citr = filter_vs_full_score_vec::const_iterator;

		} // namespace filter

		/// \brief TODOCUMENT
		using view_cache_vec = std::vector<view_cache>;

	} // namespace index

	namespace geom {

		/// \brief Type alias for a pair of coord and coord_linkage
		using coord_coord_linkage_pair         = std::pair<coord, coord_linkage>;

		/// \brief Type alias for a vector of coord_coord_linkage_pair values
		using coord_coord_linkage_pair_vec     = std::vector<coord_coord_linkage_pair>;

		/// \brief Type alias for coord_coord_linkage_pair_vec's iterator type
		using coord_coord_linkage_pair_vec_itr = coord_coord_linkage_pair_vec::iterator;

		/// \brief Type alias for a pair of two coords
		using coord_coord_pair                 = std::pair<coord, coord>;

		/// \brief Type alias for a vector of coords
		using coord_vec                        = std::vector<coord>;

		/// \brief Type alias for coord_vec's iterator type
		using coord_vec_itr                    = coord_vec::iterator;

		/// \brief TODOCUMENT
		using coord_vec_vec                    = std::vector<coord_vec>;

		/// \brief TODOCUMENT
		using coord_list_vec                   = std::vector<coord_list>;

		/// \brief TODOCUMENT
		using rotation_vec                     = std::vector<rotation>;

		/// \brief TODOCUMENT
		using coord_rot_pair                   = std::pair<coord, rotation>;

		/// \brief TODOCUMENT
		using coord_list_coord_list_pair       = std::pair<coord_list, coord_list>;

		/// \brief TODOCUMENT
		using doub_angle                       = angle<double>;

		/// \brief TODOCUMENT
		using doub_angle_doub_angle_pair       = std::pair<doub_angle, doub_angle>;

		/// \brief TODOCUMENT
		using doub_angle_vec                   = std::vector<doub_angle>;

		/// \brief TODOCUMENT
		using doub_angle_doub_angle_pair_vec   = std::vector<doub_angle_doub_angle_pair>;

		/// \brief TODOCUMENT
		template <typename T>
		using angle_angle_pair                 = std::pair<angle<T>, angle<T>>;

		/// \brief TODOCUMENT
		template <typename T>
		using angle_vec                        = std::vector<angle<T>>;

		/// \brief TODOCUMENT
		template <typename T>
		using angle_angle_pair_vec             = std::vector<angle_angle_pair<T>>;

		/// \brief TODOCUMENT
		template <typename T>
		using quat_rot_vec                     = std::vector<quat_rot_impl<T>>;
	} // namespace geom


	/// \brief TODOCUMENT
	using residue_vec                     = std::vector<residue>;

	/// \brief TODOCUMENT
	using residue_vec_itr                 = residue_vec::iterator;

	/// \brief TODOCUMENT
	using residue_vec_citr                = residue_vec::const_iterator;


	/// \brief TODOCUMENT
	using protein_vec                     = std::vector<protein>;

	/// \brief TODOCUMENT
	using sec_struc_planar_angles_vec     = std::vector<sec_struc_planar_angles>;

	/// \brief TODOCUMENT
	using sec_struc_planar_angles_vec_vec = std::vector<sec_struc_planar_angles_vec>;


	/// \brief TODOCUMENT
	using sec_struc_vec                   = std::vector<sec_struc>;

	/// \brief TODOCUMENT
	using sec_struc_vec_itr               = sec_struc_vec::iterator;

	/// \brief TODOCUMENT
	using sec_struc_vec_citr              = sec_struc_vec::const_iterator;


	/// \brief TODOCUMENT
	using amino_acid_vec                  = std::vector<amino_acid>;

	/// \brief TODOCUMENT
	using amino_acid_vec_vec              = std::vector<amino_acid_vec>;

	/// \brief TODOCUMENT
	using amino_amino_pair                = std::pair<amino_acid, amino_acid>;

	/// \brief TODOCUMENT
	using amino_amino_pair_diff_map       = std::map<amino_amino_pair, ptrdiff_t>;

	/// \brief TODOCUMENT
	using amino_diff_vec_pair             = std::pair<amino_acid, diff_vec>;

	/// \brief TODOCUMENT
	using amino_diff_vec_pair_vec         = std::vector<amino_diff_vec_pair>;

} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_STRUCTURE_TYPE_ALIASES_HPP
