/// \file
/// \brief The selection_policy_acquirer class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_SELECTION_POLICY_ACQUIRER_SELECTION_POLICY_ACQUIRER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_SELECTION_POLICY_ACQUIRER_SELECTION_POLICY_ACQUIRER_HPP

#include "cath/alignment/alignment.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/structure/structure_type_aliases.hpp"

#include <memory>
#include <utility>

// clang-format off
namespace cath::align { class common_atom_selection_policy; }
namespace cath::align { class common_residue_selection_policy; }
namespace cath::file { class pdb; }
namespace cath::file { class pdb_list; }
namespace cath::geom { class coord_list; }
namespace cath::opts { class alignment_input_spec; }
namespace cath::opts { class cath_refine_align_options; }
// clang-format on

namespace cath::opts {

	/// \brief TODOCUMENT
	class selection_policy_acquirer final {
		/// \brief TODOCUMENT
		std::shared_ptr<align::common_residue_selection_policy> comm_res_seln_pol_ptr;

		/// \brief TODOCUMENT
		std::shared_ptr<align::common_atom_selection_policy> comm_atom_seln_pol_ptr;

		[[nodiscard]] const align::common_residue_selection_policy &get_comm_res_seln_pol_ptr_cref() const;
		[[nodiscard]] const align::common_atom_selection_policy &   get_comm_atom_seln_pol_ptr_cref() const;

	  public:
		selection_policy_acquirer(const align::common_residue_selection_policy &,
		                          const align::common_atom_selection_policy &);

		[[nodiscard]] geom::coord_list_coord_list_pair get_common_coords( const align::alignment &,
		                                                                  const file::pdb &,
		                                                                  const file::pdb &,
		                                                                  // const str_veSc &,
		                                                                  const size_t &,
		                                                                  const size_t & ) const;
		[[nodiscard]] std::string                      get_descriptive_name() const;
	};

	selection_policy_acquirer get_selection_policy_acquirer(const alignment_input_spec &);
	selection_policy_acquirer get_selection_policy_acquirer(const cath_refine_align_options &);

	geom::coord_list_coord_list_pair get_common_coords(const selection_policy_acquirer &,
	                                                   const align::alignment &,
	                                                   const file::pdb_list &,
	                                                   // const str_vec &,
	                                                   const size_t &,
	                                                   const size_t &);

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_SELECTION_POLICY_ACQUIRER_SELECTION_POLICY_ACQUIRER_HPP
