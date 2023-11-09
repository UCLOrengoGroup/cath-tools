/// \file
/// \brief The align_based_superposition_acquirer class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_SUPERPOSITION_ACQUIRER_ALIGN_BASED_SUPERPOSITION_ACQUIRER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_SUPERPOSITION_ACQUIRER_ALIGN_BASED_SUPERPOSITION_ACQUIRER_HPP

#include "cath/acquirer/selection_policy_acquirer/selection_policy_acquirer.hpp"
#include "cath/acquirer/superposition_acquirer/superposition_acquirer.hpp"
#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/file/strucs_context.hpp"

// clang-format off
namespace cath::align { class alignment; }
namespace cath::file { class pdb_list; }
namespace cath::opts { class selection_policy_acquirer; }
namespace cath::sup { class superposition; }
// clang-format on

namespace cath::opts {

	/// \brief TODOCUMENT
	class align_based_superposition_acquirer final : public cath::opts::superposition_acquirer {
	private:
		/// \brief TODOCUMENT
		std::reference_wrapper<const align::alignment> the_alignment_ref;

		/// \brief TODOCUMENT
		std::reference_wrapper<const size_size_pair_vec> spanning_tree_ref;

		/// \brief TODOCUMENT
		std::reference_wrapper<const file::strucs_context> context_ref;

		/// \brief TODOCUMENT
		selection_policy_acquirer the_selection_policy_acquirer;

		sup::superposition_context do_get_superposition(std::ostream &) const final;

	public:
		align_based_superposition_acquirer(const align::alignment &,
		                                   const size_size_pair_vec &,
		                                   const file::strucs_context &,
		                                   selection_policy_acquirer);

		align_based_superposition_acquirer(const align::alignment &&,
		                                   const size_size_pair_vec &,
		                                   const file::strucs_context &,
		                                   const selection_policy_acquirer &) = delete;

		align_based_superposition_acquirer(const align::alignment &,
		                                   const size_size_pair_vec &&,
		                                   const file::strucs_context &,
		                                   const selection_policy_acquirer &) = delete;

		align_based_superposition_acquirer(const align::alignment &,
		                                   const size_size_pair_vec &,
		                                   const file::strucs_context &&,
		                                   const selection_policy_acquirer &) = delete;

		[[nodiscard]] const align::alignment &         get_alignment() const;
		[[nodiscard]] const size_size_pair_vec &       get_spanning_tree() const;
		[[nodiscard]] const file::strucs_context &     get_strucs_context() const;
		[[nodiscard]] const selection_policy_acquirer &get_selection_policy_acquirer() const;
	};

	const file::pdb_list & get_pdbs(const align_based_superposition_acquirer &);
	const file::name_set_list & get_name_sets(const align_based_superposition_acquirer &);
	const chop::region_vec_opt_vec & get_regions(const align_based_superposition_acquirer &);

	file::pdb_list get_restricted_pdbs(const align_based_superposition_acquirer &);

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_SUPERPOSITION_ACQUIRER_ALIGN_BASED_SUPERPOSITION_ACQUIRER_HPP
