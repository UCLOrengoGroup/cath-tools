/// \file
/// \brief The file_list_pdbs_acquirer class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ACQUIRER_PDBS_ACQUIRER_FILE_LIST_PDBS_ACQUIRER_H
#define _CATH_TOOLS_SOURCE_UNI_ACQUIRER_PDBS_ACQUIRER_FILE_LIST_PDBS_ACQUIRER_H

#include <boost/filesystem.hpp>

#include "acquirer/pdbs_acquirer/pdbs_acquirer.hpp"
#include "common/path_type_aliases.hpp"

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class file_list_pdbs_acquirer final : public pdbs_acquirer  {
		private:
			/// \brief TODOCUMENT
			path_vec files;

			std::unique_ptr<pdbs_acquirer> do_clone() const final;
			file::pdb_list_name_set_list_pair do_get_pdbs_and_names(std::istream &) const final;

		public:
			explicit file_list_pdbs_acquirer(path_vec);
		};

	} // namespace opts
} // namespace cath

#endif
