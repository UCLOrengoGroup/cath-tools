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

#ifndef _CATH_TOOLS_SOURCE_ACQUIRER_SUPERPOSITION_ACQUIRER_ALIGN_BASED_SUPERPOSITION_ACQUIRER_H
#define _CATH_TOOLS_SOURCE_ACQUIRER_SUPERPOSITION_ACQUIRER_ALIGN_BASED_SUPERPOSITION_ACQUIRER_H

#include <boost/filesystem.hpp>

#include "acquirer/selection_policy_acquirer/selection_policy_acquirer.h"
#include "acquirer/superposition_acquirer/superposition_acquirer.h"

namespace cath { namespace align { class alignment; } }
namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace opts { class selection_policy_acquirer; } }
namespace cath { namespace sup { class superposition; } }

namespace cath {
	namespace opts {
	
		/// \brief TODOCUMENT
		class align_based_superposition_acquirer final : public cath::opts::superposition_acquirer {
		private:
			const file::pdb_list     &pdbs;
			const str_vec            &names;
			const align::alignment   &the_alignment;
			const size_size_pair_vec &spanning_tree;

			selection_policy_acquirer the_selection_policy_acquirer;

			const file::pdb_list & get_pdbs_cref() const;
			const str_vec & get_names_cref() const;

			virtual sup::superposition_context do_get_superposition(std::ostream &) const override final;

		public:
			align_based_superposition_acquirer(const file::pdb_list &,
			                                   const str_vec &,
			                                   const align::alignment &,
			                                   const size_size_pair_vec &,
			                                   const selection_policy_acquirer &);
		};

		sup::superposition hacky_multi_ssap_fuction(const file::pdb_list &,
		                                            const str_vec &,
		                                            const size_size_pair_vec &,
		                                            const boost::filesystem::path &,
		                                            const selection_policy_acquirer &,
		                                            std::ostream &);

	} // namespace opts
} // namespace cath

#endif
