/// \file
/// \brief The common_atom_select_ca_policy class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY_COMMON_ATOM_SELECT_CA_POLICY_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY_COMMON_ATOM_SELECT_CA_POLICY_HPP

#include "cath/alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"

namespace cath {
	namespace align {

		/// \brief Concrete common_atom_selection_policy that selects the two alpha carbons in the two residues
		class common_atom_select_ca_policy final : public common_atom_selection_policy {
		private:
			void do_append_common_atoms_to_coord_lists(geom::coord_list_coord_list_pair &,
			                                           const residue &,
			                                           const residue &) const final;
			std::string do_get_descriptive_name() const final;
			std::unique_ptr<common_atom_selection_policy> do_clone() const final;

			bool do_less_than_with_same_dynamic_type(const common_atom_selection_policy &) const final;
		};

	} // namespace align
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY_COMMON_ATOM_SELECT_CA_POLICY_HPP
