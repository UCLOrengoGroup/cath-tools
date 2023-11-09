/// \file
/// \brief The superpose_orient class header

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
/// MERCHANTABILITY or ORIENTNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPERPOSE_ORIENT_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPERPOSE_ORIENT_HPP

#include "cath/structure/structure_type_aliases.hpp"

// clang-format off
namespace cath::align { class alignment; }
namespace cath::file { class pdb_list; }
namespace cath::geom { class coord_list; }
namespace cath::sup { class superposition; }
// clang-format on

namespace cath::sup {

	namespace detail {

		template <typename Fn>
		geom::coord_list get_superposed_filtered_coords_in_aln_order(const superposition &,
		                                                             const align::alignment &,
		                                                             const file::pdb_list &,
		                                                             Fn &&);

	} // namespace detail

	geom::coord_list get_coords_to_orient(const superposition &,
	                                      const align::alignment &,
	                                      const file::pdb_list &);

	geom::coord_rot_pair get_orienting_transformation(const superposition &,
	                                                  const align::alignment &,
	                                                  const file::pdb_list &);

	void orient_superposition(superposition &,
	                          const geom::coord_list &);

	superposition orient_superposition_copy(superposition,
	                                        const geom::coord_list &);


	void orient_superposition(superposition &,
	                          const align::alignment &,
	                          const file::pdb_list &);

	superposition orient_superposition_copy(superposition,
	                                        const align::alignment &,
	                                        const file::pdb_list &);

} // namespace cath::sup

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPERPOSE_ORIENT_HPP
