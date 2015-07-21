/// \file
/// \brief The pdb_input_options_block class header

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

#ifndef PDB_INPUT_OPTIONS_BLOCK_H_INCLUDED
#define PDB_INPUT_OPTIONS_BLOCK_H_INCLUDED

#include <boost/ptr_container/ptr_vector.hpp>

#include "options/options_block/options_block.h"

namespace cath { namespace opts { class pdbs_acquirer; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class pdb_input_options_block final : public cath::opts::options_block {
		private:
			static const std::string PO_PDB_INFILE;
			static const std::string PO_PDBS_FROM_STDIN;

			path_vec input_files;
			bool read_from_stdin;

			virtual std::unique_ptr<options_block> do_clone() const override final;
			virtual std::string do_get_block_name() const override final;
			virtual void do_add_visible_options_to_description(boost::program_options::options_description &) override final;
			virtual opt_str do_invalid_string() const override final;

			const path_vec & get_input_files_cref() const;
			bool get_read_from_stdin() const;

		public:
			virtual ~pdb_input_options_block() noexcept = default;

			boost::ptr_vector<pdbs_acquirer> get_pdbs_acquirers() const;
		};

	}
}

#endif
