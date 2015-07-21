/// \file
/// \brief The check_pdb_options_block class header

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

#ifndef CHECK_PDB_OPTIONS_BLOCK_H_INCLUDED
#define CHECK_PDB_OPTIONS_BLOCK_H_INCLUDED

#include "options/options_block/options_block.h"

namespace cath {
	namespace opts {

		/// \brief Handle the options specific to check-pdb
		class check_pdb_options_block final : public options_block {
		private:
			using super = options_block;

			/// \brief The PDB file to check
			boost::filesystem::path pdb_file;

			/// \brief Whether or not to permit no ATOM records
			bool permit_no_atoms;

			virtual std::unique_ptr<options_block> do_clone() const override final;
			virtual std::string do_get_block_name() const override final;
			virtual void do_add_visible_options_to_description(boost::program_options::options_description &) override final;
			virtual opt_str do_invalid_string() const override final;

		public:
			explicit check_pdb_options_block();
			virtual ~check_pdb_options_block() noexcept = default;

			boost::filesystem::path get_pdb_file() const;
			bool get_permit_no_atoms() const;

			static const std::string PO_PDB_FILE;
			static const std::string PO_PERMIT;
		};
		
	}
}

#endif
