/// \file
/// \brief The view_cache class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_VIEW_CACHE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_VIEW_CACHE_HPP

#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/structure_type_aliases.hpp"

// clang-format off
namespace cath { class protein; }
// clang-format on

namespace cath::index {

	/// \brief Cache of views (ie vectors implemented as coords) between pairs of residues in a particular list
	///        (most likely a protein)
	class view_cache final {
	private:
		/// \brief The views - coords indexed by the from-residue and then by the to-residue
		geom::coord_vec_vec views;

		static geom::coord_vec_vec build_views(const protein &);

	public:
		explicit view_cache(const protein &);

		[[nodiscard]] const geom::coord &get_view( const size_t &, const size_t & ) const;
	};

	/// \brief Getter for the view from residue with the specified from-index to the residue with the specified to-index
	inline const geom::coord & view_cache::get_view(const size_t &prm_from_index, ///< The index of the from-residue of the view to be retrieved
	                                                const size_t &prm_to_index    ///< The index of the to-residue   of the view to be retrieved
	                                                ) const {
		return views[ prm_from_index ][ prm_to_index   ];
	}

} // namespace cath::index

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_VIEW_CACHE_HPP
