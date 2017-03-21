/// \file
/// \brief The pdb_input_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_OPTIONS_OPTIONS_BLOCK_PDB_INPUT_OPTIONS_BLOCK_H
#define _CATH_TOOLS_SOURCE_OPTIONS_OPTIONS_BLOCK_PDB_INPUT_OPTIONS_BLOCK_H

#include "options/options_block/options_block.hpp"
#include "options/options_block/pdb_input_spec.hpp"

namespace cath {
	namespace opts {

		/// \brief An options block for specifying how PDBs should be read in
		class pdb_input_options_block final : public cath::opts::options_block {
		private:
			static const std::string PO_PDB_INFILE;
			static const std::string PO_PDBS_FROM_STDIN;

			/// \brief The pdb_input_spec to be configured by this options block
			pdb_input_spec the_pdb_input_spec;

			std::unique_ptr<options_block> do_clone() const final;
			std::string do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &) final;
			str_opt do_invalid_string(const boost::program_options::variables_map &) const final;

		public:
			const pdb_input_spec & get_pdb_input_spec() const;
		};

		size_t get_num_acquirers(const pdb_input_options_block &);

	} // namespace opts
} // namespace cath

#endif
