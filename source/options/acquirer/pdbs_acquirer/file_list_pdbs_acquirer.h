/// \file
/// \brief The file_list_pdbs_acquirer class header

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

#ifndef FILE_LIST_PDBS_ACQUIRER_H_INCLUDED
#define FILE_LIST_PDBS_ACQUIRER_H_INCLUDED

#include <boost/filesystem.hpp>

#include "options/acquirer/pdbs_acquirer/pdbs_acquirer.h"

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class file_list_pdbs_acquirer final : public pdbs_acquirer  {
		private:
			/// \brief TODOCUMENT
			path_vec files;

			virtual std::unique_ptr<pdbs_acquirer> do_clone() const override final;
			virtual file::pdb_list_str_vec_pair do_get_pdbs_and_names(std::istream &) const override final;

		public:
			file_list_pdbs_acquirer(const path_vec &);
			virtual ~file_list_pdbs_acquirer() noexcept = default;
		};

	}
}

#endif