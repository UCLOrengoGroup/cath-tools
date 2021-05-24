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

#ifndef _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_PDB_INPUT_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_PDB_INPUT_OPTIONS_BLOCK_HPP

#include <string_view>

#include "cath/options/options_block/options_block.hpp"
#include "cath/options/options_block/pdb_input_spec.hpp"

namespace cath::opts {

	/// \brief An options block for specifying how PDBs should be read in
	class pdb_input_options_block final : public cath::opts::options_block {
	private:
		/// \brief The pdb_input_spec to be configured by this options block
		pdb_input_spec the_pdb_input_spec;

		[[nodiscard]] std::unique_ptr<options_block> do_clone() const final;
		[[nodiscard]] std::string                    do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

	  public:
		[[nodiscard]] const pdb_input_spec &get_pdb_input_spec() const;

		/// \brief The option name for the a list of PDB files that should be read
		static constexpr ::std::string_view PO_PDB_INFILE{ "pdb-infile" };

		/// \brief The option name for whether to read PDBs from stdin
		static constexpr ::std::string_view PO_PDBS_FROM_STDIN{ "pdbs-from-stdin" };
	};

	size_t get_num_acquirers(const pdb_input_options_block &);

} // namespace cath::opts

#endif // _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_PDB_INPUT_OPTIONS_BLOCK_HPP
