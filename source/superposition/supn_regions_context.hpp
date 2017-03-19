/// \file
/// \brief The supn_regions_context class header

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

#ifndef _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPN_REGIONS_CONTEXT_H
#define _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPN_REGIONS_CONTEXT_H

#include "chopping/chopping_type_aliases.hpp"

namespace cath { namespace sup { class superposition_content_spec; } }

namespace cath {
	namespace sup {

		/// \brief The context to include in superpositions when showing the specified region(s) of structure
		enum class supn_regions_context : char {
			ALONE,    ///< Only include the specified region(s) of structure
			IN_CHAIN, ///< Include the full chain of the specified region(s) of structure
			IN_PDB    ///< Include the full PDB of the specified region(s) of structure
		};

		chop::region_vec_opt get_regions_expanded_for_context(const chop::region_vec &,
		                                                      const supn_regions_context &);

		chop::region_vec_opt get_regions_expanded_for_context(const chop::region_vec &,
		                                                      const superposition_content_spec &);

	} // namespace rslv
} // namespace cath

#endif
