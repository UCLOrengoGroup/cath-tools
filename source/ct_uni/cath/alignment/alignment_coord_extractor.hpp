/// \file
/// \brief The alignment_coord_extractor class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGNMENT_COORD_EXTRACTOR_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGNMENT_COORD_EXTRACTOR_HPP

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/common_atom_selection_policy/common_atom_select_ca_policy.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "cath/structure/structure_type_aliases.hpp"

#include <utility>
#include <vector>

namespace cath { namespace align { class common_atom_selection_policy;    } }
namespace cath { namespace align { class common_residue_selection_policy; } }
namespace cath { namespace file { class pdb; } }
//namespace cath { namespace geom { class coord_list; } }
namespace cath { class protein;    }
namespace cath { class residue;    }

namespace cath {
	namespace align {

		/// \brief Extract common coordinates from a pair of structures where common coordinates are defined by an alignment
		///
		/// \todo Refactor as much common code as possible out of these two methods
		///
		/// ATM, the first is basically a chunk of SSAP code that has been moved here
		/// and the second is a new subroutine to do similar work for residue_name_aligner alignments in cath_superpose.
		///
		/// ATM, one difference is that the first can be configured with a prm_min_sup_score whereas
		/// in the second case, the alignment will never be scored so this is of no use.
		///
		/// In the future, it might be nice to be able to score residue_name_aligner alignments
		/// and then use those scores here.
		///
		/// \todo Alter alignment to "know" whether it has been scored so that an error can be thrown
		///       on an attempt to use scores from an unscored alignment
		class alignment_coord_extractor final {
		private:
			alignment_coord_extractor() = delete;

			using residue_cref                       = std::reference_wrapper<const residue>;
			using residue_cref_residue_cref_pair     = std::pair<residue_cref, residue_cref>;
			using residue_cref_residue_cref_pair_vec = std::vector<residue_cref_residue_cref_pair>;

			static residue_cref_residue_cref_pair_vec get_common_residues(const alignment &,
			                                                              const protein &,
			                                                              const protein &,
			                                                              const common_residue_selection_policy & = common_residue_select_all_policy(),
			                                                              const alignment::size_type            & = alignment::PAIR_A_IDX,
			                                                              const alignment::size_type            & = alignment::PAIR_B_IDX);

		public:
			static std::pair<geom::coord_list_vec, geom::coord_list_vec> get_common_coords_by_residue(const alignment &,
			                                                                                          const protein &,
			                                                                                          const protein &,
			                                                                                          const common_residue_selection_policy & = common_residue_select_all_policy(),
			                                                                                          const common_atom_selection_policy    & = common_atom_select_ca_policy(),
			                                                                                          const alignment::size_type            & = alignment::PAIR_A_IDX,
			                                                                                          const alignment::size_type            & = alignment::PAIR_B_IDX);

			static geom::coord_list_coord_list_pair get_common_coords( const alignment &,
			                                                           const protein &,
			                                                           const protein &,
			                                                           const common_residue_selection_policy & = common_residue_select_all_policy(),
			                                                           const common_atom_selection_policy &    = common_atom_select_ca_policy(),
			                                                           const alignment::size_type &            = alignment::PAIR_A_IDX,
			                                                           const alignment::size_type &            = alignment::PAIR_B_IDX );

			static geom::coord_list_coord_list_pair get_common_coords( const alignment &,
			                                                           const file::pdb &,
			                                                           const file::pdb &,
			                                                           const common_residue_selection_policy & = common_residue_select_all_policy(),
			                                                           const common_atom_selection_policy &    = common_atom_select_ca_policy(),
			                                                           const alignment::size_type &            = alignment::PAIR_A_IDX,
			                                                           const alignment::size_type &            = alignment::PAIR_B_IDX,
			                                                           const ostream_ref_opt &                 = ::std::nullopt );

		};
	} // namespace align
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGNMENT_COORD_EXTRACTOR_HPP
