/// \file
/// \brief The pdbs_acquirer class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ACQUIRER_PDBS_ACQUIRER_PDBS_ACQUIRER_HPP
#define _CATH_TOOLS_SOURCE_UNI_ACQUIRER_PDBS_ACQUIRER_PDBS_ACQUIRER_HPP

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/chopping/domain/domain.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/file_type_aliases.hpp"

#include <memory>
#include <utility>

namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace file { class strucs_context; } }
namespace cath { namespace opts { class cath_refine_align_options; } }
namespace cath { namespace opts { class cath_score_align_options; } }
namespace cath { namespace opts { class cath_superpose_options; } }
namespace cath { namespace opts { class pdb_input_options_block; } }
namespace cath { namespace opts { class pdb_input_spec; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class pdbs_acquirer {
		private:
			/// \brief Pure virtual method with which each concrete pdbs_acquirer must define how to create a clone of itself
			virtual std::unique_ptr<pdbs_acquirer> do_clone() const = 0;
			
			/// \brief TODOCUMENT
			virtual file::pdb_list_name_set_list_pair do_get_pdbs_and_names(std::istream &) const = 0;

		public:
			pdbs_acquirer() = default;
			std::unique_ptr<pdbs_acquirer> clone() const;
			virtual ~pdbs_acquirer() noexcept = default;

			pdbs_acquirer(const pdbs_acquirer &) = default;
			pdbs_acquirer(pdbs_acquirer &&) noexcept = default;
			pdbs_acquirer & operator=(const pdbs_acquirer &) = default;
			pdbs_acquirer & operator=(pdbs_acquirer &&) noexcept = default;

			file::pdb_list_name_set_list_pair get_pdbs_and_names(std::istream &,
			                                                     const bool &) const;
		};

		uptr_vec<pdbs_acquirer> get_pdbs_acquirers(const pdb_input_spec &);
		uptr_vec<pdbs_acquirer> get_pdbs_acquirers(const pdb_input_options_block &);

		std::unique_ptr<pdbs_acquirer> get_pdbs_acquirer(const pdb_input_spec &);

		file::strucs_context combine_acquired_pdbs_and_names_with_ids_and_domains(file::pdb_list,
		                                                                          file::name_set_list,
		                                                                          str_vec,
		                                                                          chop::domain_vec);

		file::strucs_context get_strucs_context(const pdbs_acquirer &,
		                                        std::istream &,
		                                        const bool &,
		                                        const str_vec & = str_vec{},
		                                        const chop::domain_vec & = chop::domain_vec{});

	} // namespace opts
} // namespace cath

#endif
