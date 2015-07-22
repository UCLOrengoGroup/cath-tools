/// \file
/// \brief The common_residue_select_all_policy class header

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

#ifndef COMMON_RESIDUE_SELECT_ALL_POLICY_H_INCLUDED
#define COMMON_RESIDUE_SELECT_ALL_POLICY_H_INCLUDED

#include "alignment/common_residue_selection_policy/common_residue_selection_policy.h"

namespace cath {
	namespace align {

		/// \brief TODOCUMENT
		class common_residue_select_all_policy final : public common_residue_selection_policy {
		private:
			virtual size_vec do_select_common_residues(const alignment &,
			                                           const std::vector<alignment::size_type> &,
			                                           const alignment::size_type &,
			                                           const alignment::size_type &) const override final;
			virtual std::string do_get_descriptive_name() const override final;
			virtual std::unique_ptr<common_residue_selection_policy> do_clone() const override final;

			virtual bool do_less_than_with_same_dynamic_type(const common_residue_selection_policy &) const override final;

		public:
			virtual ~common_residue_select_all_policy() noexcept = default;
		};
	}
}

#endif